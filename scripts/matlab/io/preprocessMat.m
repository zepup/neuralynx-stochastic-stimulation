function collData = preprocessMat(pt, sessionDate, timeStart, timeStop)
%PREPROCESSMAT Load the bundled Neuralynx example into one MATLAB structure.
% Main purpose:
% - read the small example folder under examples/neuralynx_test/1 config loop
% - load simple `.ncs` channel files and `Events.nev`
% - keep the real patient-data path/anatomy workflow as commented guidance
%
% Example-data usage:
%   addpath('/path/to/repository/scripts/matlab/io')
%   collData = preprocessMat();
%
% Optional folder override with the same simple layout:
%   setenv('NCS_PROJECT_DATA_DIR', '/path/to/folder/with/ncs/nev/cfg')
%   collData = preprocessMat('PATIENT_ID', 'YYYY-MM-DD', 'HH:MM:SS', 'HH:MM:SS');

if nargin < 1 || isempty(pt)
    pt = 'example_neuralynx_test';
end

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'shared'))
paths = project_paths();
if ~isempty(paths.fieldtripDir)
    addpath(paths.fieldtripDir)
    ft_defaults
end

dataPath = paths.dataDir;
saveDir = paths.outputDir;

%% Optional real-data anatomy workflow
% The bundled example has only two `.ncs` files, one `Events.nev`, and one
% `Condition_1.cfg`, so it does not need an anatomy spreadsheet.
%
% For real patient data, uncomment/adapt the block below after setting
% NCS_PROJECT_ANATOMY_FILE or placing `*_2mm.xlsx` inside the data folder.
%
% containsRemove = {'.', '..', 'csf', 'choroid-plexus', ...
%     'hypointensities', 'unknown', 'Lateral-Ventricle'};
% electrodeTable = readtable(paths.anatomyFile);
% contactList = electrodeTable(~contains(electrodeTable.Unit, containsRemove), :);

%% Find example Neuralynx files
ncsFiles = dir(fullfile(dataPath, '*.ncs'));
if isempty(ncsFiles)
    error('No .ncs files found in %s', dataPath);
end

channelNames = erase({ncsFiles.name}, '.ncs')';
contactList = table(channelNames, 'VariableNames', {'contact'});

eventFiles = dir(fullfile(dataPath, 'Events*.nev'));
eventsColl = table();
if ~isempty(eventFiles)
    eventName = fullfile(dataPath, eventFiles(1).name);
    eventsColl = struct2table(ft_read_event(eventName));
end

cfgFiles = dir(fullfile(dataPath, 'Condition_*.cfg'));
cfgText = "";
if ~isempty(cfgFiles)
    cfgText = string(fileread(fullfile(dataPath, cfgFiles(1).name)));
end

%% Optional real-data time crop
cropToWindow = nargin >= 4 && ~isempty(sessionDate) && ~isempty(timeStart) && ~isempty(timeStop);
if cropToWindow
    unixStart = uint64(1e6 * posixtime(datetime([sessionDate ' ' timeStart])));
    unixEnd = uint64(1e6 * posixtime(datetime([sessionDate ' ' timeStop])));
    if ~isempty(eventsColl) && ismember('timestamp', eventsColl.Properties.VariableNames)
        eventsColl = eventsColl(eventsColl.timestamp >= unixStart & eventsColl.timestamp <= unixEnd, :);
    end
else
    unixStart = uint64(0);
    unixEnd = uint64(intmax('uint64'));
end

%% Load channel data
globalData = cell(numel(ncsFiles), 1);
globalTime = cell(numel(ncsFiles), 1);

for i = 1:numel(ncsFiles)
    loadName = fullfile(dataPath, ncsFiles(i).name);
    data = ft_read_data(loadName);
    hdr = ft_read_header(loadName);

    data = data * hdr.orig.ADBitVolts;
    curTime = (((1:length(data)) - 1) / hdr.Fs);
    curTime = hdr.orig.FirstTimeStamp + uint64(curTime * 1e6);

    keepIdx = curTime >= unixStart & curTime <= unixEnd;
    globalTime{i} = curTime(keepIdx);
    globalData{i} = data(keepIdx);
end

%% Compile data
collData = struct();
collData.dataDir = dataPath;
collData.time = globalTime;
collData.events = eventsColl;
collData.contactList = contactList;
collData.data = globalData;
collData.conditionConfig = cfgText;

save(fullfile(saveDir, [pt '.mat']), 'collData', '-v7.3')

if exist('struct2h5', 'file') == 2
    struct2h5(fullfile(saveDir, [pt '.h5']), collData);
else
    warning('struct2h5 is not on the MATLAB path; skipped HDF5 export.');
end
end
