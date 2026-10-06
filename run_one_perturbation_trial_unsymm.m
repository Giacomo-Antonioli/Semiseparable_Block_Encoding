

function [err,bound_val,bound_ok,u,v,x,y,u_tilde,v_tilde,x_tilde,y_tilde, ...
    S_exact,S_tilde,w_lo,w_up,alpha] = ...
    run_one_perturbation_trial_unsymm(N, epsilon)
% RUN_ONE_PERTURBATION_TRIAL_UNSYMM
%   Generates a random unsymmetric semiseparable problem, perturbs the
%   (normalized) generators u,v,x,y by vectors of norm epsilon, rebuilds
%   the block-encoding circuit with the perturbed generators, extracts
%   the encoded matrix S_tilde by simulation, and compares it against
%   the exact (unperturbed) target matrix S_exact.
%
%   The reference bound follows Theorem thm:unsym:
%
%       w_lo  = nu*nv,  w_up = nx*ny
%       alpha = sqrt(N)*(w_lo + w_up)
%       || S - alpha*S_tilde_raw ||_2 <= 2*sqrt(w_lo^2+w_up^2)*epsilon
%
%   Since S_exact = S and S_tilde IS the raw circuit output
%   (S_tilde_raw), the bound in this (already normalized) comparison
%   scale becomes:
%
%       || S_exact - S_tilde*alpha ||_2 <= (2*sqrt(w_lo^2+w_up^2))*epsilon

    % --- exact (unperturbed) problem -----------------------------
    [S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N);

    u = u/nu;
    v = v/nv;
    x = x/nx;
    y = y/ny;

    w_lo  = nu*nv;
    w_up  = nx*ny;
    alpha = sqrt(N)*(w_lo+w_up);

    S_exact = S;

    % --- perturb the (normalized) generator directions -----------
    delta_u = randn(N,1); delta_u = epsilon*delta_u/norm(delta_u);
    delta_v = randn(N,1); delta_v = epsilon*delta_v/norm(delta_v);
    delta_x = randn(N,1); delta_x = epsilon*delta_x/norm(delta_x);
    delta_y = randn(N,1); delta_y = epsilon*delta_y/norm(delta_y);

    u_tilde = u + delta_u;
    v_tilde = v + delta_v;
    x_tilde = x + delta_x;
    y_tilde = y + delta_y;

    u_tilde =u_tilde/norm(u_tilde);
    v_tilde =v_tilde/norm(v_tilde);
    x_tilde =x_tilde/norm(x_tilde);
    y_tilde =y_tilde/norm(y_tilde);
    % --- build circuit with the perturbed generators --------------
    [circ,~,~] = build_unsymmetric_semiseparable_circuit( ...
        u_tilde, v_tilde, x_tilde, y_tilde, nu, nv, nx, ny);

    dim = 2^circ.nbQubits;
    psi = zeros(dim,1);
    S_tilde = zeros(N);

    for j = 1:N
        psi(:) = 0;
        psi(j) = 1;
        sim = circ.simulate(psi);
        out = sim.states;
        S_tilde(:,j) = out(1:N);
    end

    % --- error against the exact (unperturbed) target -------------
    err = norm(S_exact - S_tilde*alpha, 2);

    % --- theorem bound (normalized scale), computed per trial -----
    bound_val = (2*sqrt(w_lo^2+w_up^2))*epsilon;
    bound_ok  = err <= bound_val;

end