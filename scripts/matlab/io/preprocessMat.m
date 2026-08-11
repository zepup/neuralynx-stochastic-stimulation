function collData = preprocessMat()
%PREPROCESSMAT Load Neuralynx `.ncs` and `.nev` files into one MATLAB structure.
% Main purpose:
% - identify recording blocks that cover a requested time range
% - load channel data and event timestamps from Neuralynx files
% - save a simple MATLAB structure for later inspection

clc; clear;

%% meta Data
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'shared'))
paths = project_paths();
if ~isempty(paths.fieldtripDir)
    addpath(paths.fieldtripDir)
    ft_defaults
end

containsRemove=[{'.'},{'..'},{'csf'}, {'choroid-plexus'}, {'hypointensities'}, {'unknown'}, {'Lateral-Ventricle'}];

%% Specific meta data 
%{
Because I have task data that has a specific time/date based on the
epilepsy tracker sheet, i added this here. NLX operates on unix time so i
needed to be able to convert those task start/end times to unix so things
were better parsable. Below contains probably a lot of pathways you will
need to change for your own directories

Notes for reuse:
- real recording suffixes and event files are often split because each
  Neuralynx recording block is commonly capped at about 4 hours
- one stimulation session can span multiple recording labels such as 0029
  and 0030
%}

pt='SUBJECT_ID';
date='YYYY-MM-DD';
timeStart='HH:MM:SS';
timeStop='HH:MM:SS';
% the following 2 lines are in unix microsecond time (same as nlx)
unixStart = uint64(1e6*posixtime(datetime([date ' ' timeStart])));
unixEnd =  uint64(1e6*posixtime(datetime([date ' ' timeStop])));
dataPath = paths.dataDir;
electrodeTable = readtable(fullfile(paths.dataDir, 'anatomy.xlsx'));
contactList = electrodeTable(~contains(electrodeTable.Unit, containsRemove),:);
saveDir = paths.outputDir;

%{
 This looks for files in the range of your task/neural data of interest.
 Not sure how it works when the task spans 2 files because i havent had to
 run into this yet. Keep that in mind. 
%}
[startLabel, endLabel]=pullTimeLabels(date, dataPath, timeStart, timeStop, contactList);
labels = unique({startLabel, endLabel});

%% Check for valid contacts
%{
Not all contacts are in brain areas we care about. This loops through the
imaging sheet and checks what files exist as well as if its in an
interesting area so you dont spend time loading data you dont need. 
%}
disp('Checking for existing .ncs files...')

validIdx = false(height(contactList), 1);
for i = 1:height(contactList)
    curContact = contactList.contact{i};
    % Check if *any* file for that contact exists for any label
    existsFlag = false;
    for j = 1:length(labels)
        loadName = fullfile(dataPath, [curContact '_' labels{j} '.ncs']);
        if isfile(loadName)
            existsFlag = true;
            break; % no need to check further labels
        end
    end
    validIdx(i) = existsFlag;
end
contactList = contactList(validIdx, :);

%% Compile events
%{
a lot of the code that was previously in here was specific to my stop
signal task. I took that out as it would create unnecessary confusion but
you would have to edit this for your specific task/stim events. 
%}
eventsColl=[];
for j = 1:length(labels)
    eventName =[dataPath 'Events_' labels{j} '.nev'];
    events = struct2table(ft_read_event(eventName));
    events=events((events.timestamp>=unixStart & events.timestamp<=unixEnd),:);
    eventsColl=[eventsColl; events];
    eventsColl=eventsColl((eventsColl.timestamp>=unixStart & eventsColl.timestamp<=unixEnd),:);
end

%% Load data using parallel processing to speed up workflow
%{
parallel processing does not necessarily need to be implemented but this is
pretty agonizingly slow. It seems parallelization helps a bit but again
still slow because matlab
%}
delete(gcp('nocreate'))
numCores = feature('numcores'); 
c = parcluster('local');           
disp(c.NumWorkers)                
parpool('local', numCores - 1)

parfor i = 1:size(contactList,1) %would change this to "for" rather than "parfor" if you want to troubleshoot it
    disp(['contact ' num2str(i) ' of ' num2str(size(contactList,1))])
    curContact = contactList.contact{i};
    collapData=[];
    collapTime=[];
    for j = 1:length(labels)
        loadName =[dataPath curContact '_' labels{j} '.ncs'];

        data=ft_read_data(loadName); %ft toolbox for reading data and header
        hdr = ft_read_header(loadName);
        data=data*hdr.orig.ADBitVolts;
        curTime=(((1:length(data))-1)/hdr.Fs);
        curTime = hdr.orig.FirstTimeStamp+uint64(curTime*1e6);

        collapData=[collapData, data];
        collapTime = [collapTime, curTime];
    end
    
    collapTime=collapTime(collapTime>=unixStart & collapTime<=unixEnd);
    collapData=collapData(collapTime>=unixStart & collapTime<=unixEnd);

    globalTime(i,:)=collapTime;
    globalData(i,:)=collapData;
end

%% Compile data
collData=struct();
collData.time=globalTime(1,:);
collData.events=eventsColl;
collData.contactList=(contactList);
collData.data=globalData;
save([saveDir pt '.mat'], 'collData', '-v7.3')

%%
struct2h5([saveDir pt '.h5'], collData);
end
