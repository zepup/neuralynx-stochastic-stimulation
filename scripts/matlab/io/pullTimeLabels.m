function [startLabel, endLabel]=pullTimeLabels(date, dataPath, timeStart, timeEnd, contactList)
%PULLTIMELABELS Find the first and last Neuralynx recording suffixes that
%cover a requested time interval.
% Main purpose:
% - inspect available `.ncs` files for one contact
% - map a clock-time window onto recording labels such as `0029` or `0030`

recordingNames = {dir(dataPath).name};
contactUnique = unique(recordingNames(find(contains(recordingNames,contactList.contact{1}))));

for i = 1:length(contactUnique)
    curContact = contactUnique{i};
    fullFileName = fullfile(dataPath, curContact);
    hdr = ft_read_header(fullFileName);

    t1 = datetime(hdr.orig.FirstTimeStamp/1e6, 'ConvertFrom', 'posixtime');

    if (t1<datetime([date ' ' timeStart]))==1
        startLabel=split(curContact,'_');
        startLabel=split(startLabel{2},'.ncs');
        startLabel=startLabel{1};
        break;
    end
end

for i = 1:length(contactUnique)
    curContact = contactUnique{i};
    fullFileName = fullfile(dataPath, curContact);
    hdr = ft_read_header(fullFileName);

    t2 = datetime(hdr.orig.LastTimeStamp/1e6, 'ConvertFrom', 'posixtime');

    if (t2>datetime([date ' ' timeEnd]))==1
        endLabel=split(curContact,'_');
        endLabel=split(endLabel{2},'.ncs');
        endLabel=endLabel{1};
        break;
    end
end

end
