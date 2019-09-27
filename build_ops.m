function [ G, D, I ] = build_ops( Grid )
% author: Daniel Espinoza
% date: 17-Sept-19
% description:
% This function computes the discrete divergence and gradient matrices on a
% regular staggered grid using central difference approximations. The
% discrete gradient assumes homogeneous boundary conditions.
e = ones(Grid.Nx,1);
D = spdiags( [-e e],0:1,Grid.Nx,Grid.Nfx);
G=-D';
I= spdiags(e,Grid.Nx,Grid.Nx);
G(1,1)=0;G(Grid.Nfx,Grid.Nx)=0;
end

