function paths = project_paths(dataDir)
%PROJECT_PATHS Centralized path configuration for the minimal MATLAB scripts.
% Main purpose:
% - point MATLAB to the included example Neuralynx data by default
% - let users optionally override paths for real data outside the repository
%
% Usage:
%   paths = project_paths();
%
% Optional real-data override:
%   setenv('NCS_PROJECT_DATA_DIR', '/path/to/real/session/stochasticStim')
%   setenv('NCS_PROJECT_ANATOMY_FILE', '/path/to/PATIENT_ID_2mm.xlsx')
%   paths = project_paths();

    if nargin < 1 || isempty(dataDir)
        dataDir = getenv('NCS_PROJECT_DATA_DIR');
    end

    if isempty(dataDir)
        sharedDir = fileparts(mfilename('fullpath'));
        repoRoot = fullfile(sharedDir, '..', '..', '..');
        dataDir = fullfile(repoRoot, 'examples', 'neuralynx_test', '1 config loop');
    end

    paths.dataDir = dataDir;
    if isempty(paths.dataDir)
        error(['Set NCS_PROJECT_DATA_DIR to the folder containing Neuralynx data, ' ...
               'event files, and derived metadata.']);
    end

    paths.fieldtripDir = getenv('FIELDTRIP_DIR');
    if isempty(paths.fieldtripDir)
        paths.fieldtripDir = '';
    end

    paths.outputDir = getenv('NCS_PROJECT_OUTPUT_DIR');
    if isempty(paths.outputDir)
        paths.outputDir = paths.dataDir;
    end

    % Real patient datasets may include an anatomy spreadsheet such as
    % PATIENT_ID_2mm.xlsx. The bundled example intentionally does not rely on
    % anatomy metadata; scripts should tolerate an empty anatomyFile.
    paths.anatomyFile = getenv('NCS_PROJECT_ANATOMY_FILE');
    if isempty(paths.anatomyFile)
        anatomyCandidates = dir(fullfile(paths.dataDir, '*_2mm.xlsx'));
        anatomyCandidates = anatomyCandidates(~startsWith({anatomyCandidates.name}, '~$'));
        if isempty(anatomyCandidates)
            anatomyCandidates = dir(fullfile(paths.dataDir, 'anatomy.xlsx'));
            anatomyCandidates = anatomyCandidates(~startsWith({anatomyCandidates.name}, '~$'));
        end
        if isempty(anatomyCandidates)
            paths.anatomyFile = '';
        else
            paths.anatomyFile = fullfile(paths.dataDir, anatomyCandidates(1).name);
        end
    end
end
