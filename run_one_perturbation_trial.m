    
    function [err,bound_ok,u,v,u_tilde,v_tilde,S,S_tilde] =run_one_perturbation_trial(n,type,epsilon,bound_const)
    
    N = 2^n;
    [u,v] = local_generate_uv(n,type);

    u = u(:);
    v = v(:);

    u = u/norm(u);
    v = v/norm(v);

    % random perturbations
    delta_u = randn(N,1);
    delta_u = epsilon*delta_u/norm(delta_u);

    delta_v = randn(N,1);
    delta_v = epsilon*delta_v/norm(delta_v);

    u_tilde = u + delta_u;
    v_tilde = v + delta_v;
    u_tilde=u_tilde/norm(u_tilde);
    v_tilde=v_tilde/norm(v_tilde);
    S = tril(u*v') + triu(v*u',1);
    
    circ = build_semiseparable_circuit(u_tilde,v_tilde);

    dim = 2^circ.nbQubits;
    S_tilde = zeros(N);
    
    psi = zeros(dim,1);
    
    for j = 1:N
    
        psi(:) = 0;
        psi(j) = 1;
    
        sim = circ.simulate(psi);
        out = sim.states;
    
        S_tilde(:,j) = 2*sqrt(N) * out(1:N);
    
    end

    
    err = norm(S-S_tilde,2);
    bound_ok = err < bound_const*epsilon;
        
    end
    
    
    

