function [tot_thetas, tot_phis, tot_gammas] = get_all_angles(state)
    %GET_ALL_ANGLES Iteratively build a list of rotation angles for a state
    %
    %   [tot_thetas, tot_phis, tot_gammas] = get_all_angles(state)
    %
    %   Inputs:
    %       state - column vector of coefficients of the quantum state
    %
    %   Outputs:
    %       tot_thetas, tot_phis, tot_gammas - cell arrays of rotation angles
    
    n = log2(length(state));
    
    % Initialize cell arrays
    tot_thetas = {};
    tot_phis = {};
    tot_gammas = {};
    
    % First call to get_angles
    [thetas, phis, gammas, newstate] = get_angles(state);
    tot_thetas{end+1} = thetas;
    tot_phis{end+1} = phis;
    tot_gammas{end+1} = gammas;
    
    % Iteratively get angles for smaller states
    for j1 = 1:(n-1)
        [thetas, phis, gammas, newstate] = get_angles(newstate);
        tot_thetas{end+1} = thetas;
        tot_phis{end+1} = phis;
        tot_gammas{end+1} = gammas;
    end
    
    % Reverse the order to match iterative building from |ψ>_1 to |ψ>_n
    tot_thetas = flip(tot_thetas);
    tot_phis = flip(tot_phis);
    tot_gammas = flip(tot_gammas);
end
