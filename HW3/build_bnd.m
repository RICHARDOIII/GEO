function [B,N,fn] = build_bnd(Param,Grid,I)
%author: Daniel Espinoza
% date: 9/27/2019
% Description:% This function computes 
%the operators and r.h.s vectors for both Dirichlet% and Neumann boundary conditions.
dof_dor=[Param.dof_dir,Param.dof_neu];
B=I(dof_dor,:);
N=I;
N(:,dof_dor)=[];
fn=spalloc(Grid.Nx,1,length(Param.qb));
fn(Param.dof_neu)=Param.qb;
end

