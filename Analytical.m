function [ h ] = Analytical( Param,Grid,D,G)
%ANALYTICAL Summary of this function goes here
%   Detailed explanation goes here
    e = ones(Grid.Nfx,1);
    k=Param.K;
    Kd=spdiags(e*(1/k),0,Grid.Nfx,Grid.Nfx);
    Kd(1,1)=0;Kd(Grid.Nfx,Grid.Nfx)=0;
    
    fs=spalloc(Grid.Nx,1,1)-(Param.qp/Param.b);
    h=fs\(-D*Kd*G);
end

