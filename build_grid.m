function [ Grid ] = build_grid(Grid)
% author: Daniel Espinoza
% date: 03-Sept-2019
%BUILD_GRID :This function computes takes in minimal definition of the computational
% domain and grid and computes all containing all pertinent information
% about the grid.
% Scalar length of the domain
Grid.Lx = Grid.xmax- Grid.xmin;
% Scalar cell width
Grid.dx = Grid.Lx/Grid.Nx;
% number of fluxes
Grid.Nfx = Grid.Nx+1;
% Nx by 1 column vector of cell center locations
Grid.xc = Grid.xmin+Grid.dx/2:Grid.dx:Grid.xmax;
% Nfx by 1 column vector of cell face locations
Grid.xf = Grid.xmin:Grid.dx:Grid.xmax;
% Grid.dof = Nx by 1 column vector from 1 to N containing the degrees of
% freedom, i.e. cell numbers 
Grid.dof = 1:Grid.Nx;
% Grid.dof_f = Nx by 1 column vector from 1 to Nfx containing the degrees
% of freedom, i.e. cell faces
Grid.dof_f = 1:Grid.Nfx;
% Grid.dof_xmin = scalar cell degree of freedom corrsponding to the left
% boundary 
Grid.dof_xmin =Grid.xc(1);
% Grid.dof_xmax = scalar cell degree of freedom corrsponding to the right
% boundary 
Grid.dof_xmax =Grid.xc(Grid.Nx);
% Grid.dof_f_xmin = scalar face degree of frxceedom corrsponding to the left
% boundary
Grid.dof_f_xmin =Grid.xf(1);
% Grid.dof_f_xmax = scalar face degree of freedom corrsponding to the right
% boundary 
Grid.dof_f_xmax =Grid.xf(Grid.Nfx);
end

