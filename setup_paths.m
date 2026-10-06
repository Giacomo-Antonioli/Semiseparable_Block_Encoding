function root = setup_paths()
%SETUP_PATHS Add the repository folders and FABLE to the MATLAB path.
%
% Run once from the repository root before using any script:
%
%   >> setup_paths
%
% Output:
%   root - absolute path of the repository root.

    root = fileparts(mfilename('fullpath'));

    addpath(root);
    addpath(genpath(fullfile(root, 'src')));
    addpath(fullfile(root, 'experiments', 'utils'));
    addpath(fullfile(root, 'plotting'));

    fableDir = fullfile(root, 'external', 'fable', 'fable-qclab');
    if exist(fableDir, 'dir')
        addpath(fableDir);
    else
        warning('setup_paths:noFable', ...
            ['FABLE not found in %s. Run "git submodule update --init" ' ...
             'to use the FABLE comparison experiments.'], fableDir);
    end

    if isempty(which('qclab.QCircuit'))
        warning('setup_paths:noQCLAB', ...
            'QCLAB is not on the MATLAB path. See README.md.');
    end

end
