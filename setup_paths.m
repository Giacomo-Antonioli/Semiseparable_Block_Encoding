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

    % Prefer the bundled QCLAB fork (external/qclab, a git submodule): it has
    % the features this code needs on top of QCLAB v1.2.0.
    qclabDir = fullfile(root, 'external', 'qclab');
    if exist(fullfile(qclabDir, '+qclab'), 'dir')
        addpath(qclabDir);
    end

    if isempty(which('qclab.QCircuit'))
        warning('setup_paths:noQCLAB', ...
            ['QCLAB is not on the MATLAB path. Run "git submodule update ' ...
             '--init" or see README.md.']);
    elseif isempty(which('qclab.qgates.CSWAP'))
        warning('setup_paths:oldQCLAB', ...
            ['The QCLAB on the path lacks qgates.CSWAP (needs the fork in ' ...
             'external/qclab). Run "git submodule update --init".']);
    end

end
