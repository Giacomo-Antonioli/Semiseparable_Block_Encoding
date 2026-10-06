function U = generate_sp_block(psi)

    dim = numel(psi);
    psi = psi(:) / norm(psi);

    basis    = [psi, eye(dim, dim - 1)];
    [U, ~]   = qr(basis);

    if norm(U(:,1) - psi) > norm(-U(:,1) - psi)
        U(:,1) = -U(:,1);
    end

end