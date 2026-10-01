function plotExampleData()
%PLOTEXAMPLEDATA Plot the bundled Neuralynx example recording.
% Main purpose:
% - call preprocessMat() to load the included example files
% - plot each example `.ncs` channel in microvolts
% - mark timestamps from `Events.nev` when event records are present
% - save a PNG under examples/plots/
%
% Usage:
%   addpath('/path/to/repository/scripts/matlab/io')
%   addpath('/path/to/repository/scripts/matlab/inspection')
%   setenv('FIELDTRIP_DIR', '/path/to/fieldtrip')
%   plotExampleData()

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fullfile(scriptDir, '..', '..', '..');
addpath(fullfile(scriptDir, '..', 'io'))

collData = preprocessMat();

if isempty(collData.time) || isempty(collData.data)
    error('No channel data found in collData.');
end

firstTimestamp = min(cellfun(@(x) x(1), collData.time));
channelCount = numel(collData.data);

fig = figure('Color', 'w', 'Position', [100 100 1200 300 * channelCount]);

for channelIndex = 1:channelCount
    subplot(channelCount, 1, channelIndex)
    timeSec = (double(collData.time{channelIndex}) - double(firstTimestamp)) / 1e6;
    signalUv = collData.data{channelIndex} * 1e6;

    plot(timeSec, signalUv, 'k', 'LineWidth', 0.7)
    hold on

    if ~isempty(collData.events) && ismember('timestamp', collData.events.Properties.VariableNames)
        eventSec = (double(collData.events.timestamp) - double(firstTimestamp)) / 1e6;
        eventSec = eventSec(eventSec >= min(timeSec) & eventSec <= max(timeSec));
        for eventIndex = 1:numel(eventSec)
            yLimits = ylim;
            line([eventSec(eventIndex) eventSec(eventIndex)], yLimits, ...
                'Color', [0.8 0 0], 'LineStyle', '--', 'LineWidth', 0.8);
        end
    end

    if istable(collData.contactList) && ismember('contact', collData.contactList.Properties.VariableNames)
        channelName = string(collData.contactList.contact(channelIndex));
    else
        channelName = "Channel " + channelIndex;
    end

    title(channelName, 'Interpreter', 'none')
    ylabel('uV')
    if channelIndex == channelCount
        xlabel('Time from recording start (s)')
    end
end

if exist('sgtitle', 'file') == 2
    sgtitle('Bundled Neuralynx Example: NCS traces with NEV event markers')
end

outputDir = fullfile(repoRoot, 'examples', 'plots');
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end
outputPath = fullfile(outputDir, 'matlab_example_trace.png');
saveas(fig, outputPath)
fprintf('Wrote %s\n', outputPath);
end
