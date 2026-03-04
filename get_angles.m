function [thetas, phis, gammas, newstate] = get_angles(state_0)
    % Initialize outputs
    thetas = [];
    phis = [];
    gammas = [];
    newstate = [];
    
    n = length(state_0)/2;
    
    for j1 = 1:n
        idx = 2*j1 - 1;  % MATLAB indexing
        
        % Compute global phase
        globalphase = conj(state_0(idx)) / norm(state_0(idx));
        
        % Compute theta
        norm_pair = sqrt(norm(state_0(idx))^2 + norm(state_0(idx+1))^2);
        theta = 2 * acos(globalphase * state_0(idx) / norm_pair);
        
        % Compute gamma and phi
        gamma = log(globalphase) / 1i;
        phi = log((globalphase * state_0(idx+1)) / sin(theta/2)) / 1i;
        
        % Store real parts if imaginary part is small
        if imag(theta) <= 1e-5
            thetas(end+1) = real(theta);
        end
        
        phis(end+1) = real(phi);
        
        if imag(gamma) <= 1e-5
            gammas(end+1) = real(2*gamma);
        end
        
        newstate(end+1) = norm_pair;
    end
end
