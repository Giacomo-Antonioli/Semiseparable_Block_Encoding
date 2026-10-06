addpath("fable\fable-qclab\")
% circ=fable([1,0;0,1],'cutoff',0,false)
% circ.draw
clear;
clc;
n=2:3;
types=1:3;
loading=1:2;

models=2;
for i=1:size(n,2)
    dim=n(i);
    N=2^dim;
for j=1:size(types,2)
    type=j;
    [u,v]=local_generate_uv(dim,type);

    if type==2
        models=3;
    else
        models=2;
for k=0:models
    
if k == 0

    circ = build_semiseparable_circuit(u, v);

elseif k == 1

    circ = build_semiseparable_circuit2(u, v);

elseif k == 2

    [circ,v,u] = build_semiseparable_circuit3(u, v);

else

    error('Invalid value of k. Expected k = 1, 2, or 3.');

end

     S      = tril(u*v') + triu(v*u', 1);
     fableCirc=fable(S,"cutoff",0,false);
     my_S = zeros(N);
     Fable_S = zeros(N);
     
    
    psiMy = zeros(2^circ.nbQubits,1);
    psiFable = zeros(2^fableCirc.nbQubits,1);
    
    for col = 1:N
    
        psiMy(:) = 0;
        psiMy(col) = 1;
        
        psiFable(:) = 0;
        psiFable(col) = 1;
    
        sim = circ.simulate(psiMy);
        out = sim.states;
        my_S(:,col) = 2*sqrt(N) * out(1:N);

        simFable = fableCirc.simulate(psiFable);
        out = simFable.states;
        Fable_S(:,col) =N * out(1:N);

    
    end

end

end
end

end