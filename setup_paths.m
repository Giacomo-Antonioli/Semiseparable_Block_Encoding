function root = setup_paths()
%SETUP_PATHS Add the library and the bundled QCLAB fork to the MATLAB path.
%
% Run once before using the circuit builders:
%
%   >> setup_paths
%
% Output:
%   root - absolute path of the repository root.

    root = fileparts(mfilename('fullpath'));

    addpath(genpath(fullfile(root, 'circuits')));
    addpath(genpath(fullfile(root, 'state_prep')));

    % Bundled QCLAB fork (git submodule): QCLAB v1.2.0 plus the features this
    % code needs (CSWAP, InitializeStateVector, circuit statistics).
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
