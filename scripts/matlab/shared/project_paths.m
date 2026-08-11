function paths = project_paths()
%PROJECT_PATHS Centralized path configuration for the minimal MATLAB scripts.
% Main purpose:
% - let users define the Neuralynx data folder once
% - optionally point MATLAB to FieldTrip and an output folder

    paths.dataDir = getenv('NCS_PROJECT_DATA_DIR');
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
end
