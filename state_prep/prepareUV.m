function [circU, circV,u,v] = prepareUV(n)
% PREPAREUV  Preparazione esatta O(n)-gate, depth O(1), di
%   u_i = exp(theta*i/N),  v_i = exp(-theta*i/N),  i=0,...,N-1, N=2^n
% Entrambi sono stati prodotto (sequenze geometriche), quindi bastano
% n rotazioni RY per registro, senza gate controllati.
    theta=3.0;
    N = 2^n;
    circU = qclab.QCircuit(n,n+1);
    circV = qclab.QCircuit(n);

    for k = 0:n-1
        w       = 2^(n-1-k);              % peso del qubit k (MSB-first)
        theta_u = 2*atan( exp(theta*w/N) );
        theta_v = pi - theta_u;           % identita' arctan(x)+arctan(1/x)=pi/2

        circU.push_back( qclab.qgates.RotationY(k, theta_u) );
        circV.push_back( qclab.qgates.RotationY(k, theta_v) );
    end

    u=circU.matrix;
    u=u(:,1);
    v=circV.matrix;
    v=v(:,1);
end

