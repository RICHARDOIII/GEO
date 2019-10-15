function [q ,r] = comp_flux(D,Kd,G,h,fs,Grid,Param)
% author: Daniel Espinoza
% date:  02/10/2019
% Description:
% Computes the mass conservative fluxes across all boundaries from the
% residual of the compatability condition over the boundary cells.
% Note: Current implmentation works for all cases where one face
%       is assigend to each bnd cell. So corner cells must have
%       natual BC’s on all but one face.%

%creating Laplacian
L=-D*Kd*G;
q=-Kd*G*h;
%solving for residual
r = L*h - fs;
dof_bnd=[Param.dof_dir;Param.dof_neu];
dof_f_bnd = [Param.dof_f_dir;Param.dof_f_neu];
Lia = ismember(dof_bnd,Grid.dof_xmin);
Lix = -ismember(dof_bnd,Grid.dof_xmax);
sign=Lia+Lix;
q(dof_f_bnd)=sign.*r(dof_bnd).*Grid.V(dof_bnd)./Grid.A(dof_f_bnd);
end

