% file name: dee467_VarCo1D.m
% author: Daniel Espinoza
% date: July 20, 1923
clc, close all, clear all
%%  Analytic solution
K1=5e-3;K2=5e-5;
h1 = @(x)-.00962*x+120.0;    % analytic vol. flux
h2 = @(x)-.962*x+139.038;    % analytic vol. flux
h3 = @(x)-.00962*x+100.962;    % analytic vol. flux
q1 = @(x).00962*K1;    % analytic vol. flux
q2 = @(x).962*K2;    % analytic vol. flux
q3 = @(x).00962*K1;    % analytic head%% Finite volume solution
Grid.xmin = 0; 
Grid.xmax = 100; 
Grid.Nx = 5;
Grid = build_grid(Grid);
Grid.N=Grid.Nx; %1D
%% Create the K vector
P1=20/Grid.dx;
P2=40/Grid.dx;
K=zeros(1,Grid.N);
K(1:P1)=K1;
K(P1+1:P2)=K2;
K(P2+1:end)=K1;
%% Get respective K matrix based on Arithmetic or Harmonic mean
% p = power of the generalized mean
% 1 (arithmetic mean)
% -1 (harmonic mean)
KdA=comp_mean(K,1,Grid);
KdH=comp_mean(K,-1,Grid);
%% Boundary conditions
Param.dof_dir =Grid.dof_xmin; 
Param.dof_f_dir =Grid.dof_f_xmin; 
% cell and face on Dirichlet bnd
Param.dof_neu=Grid.dof_xmax;
Param.dof_f_neu=Grid.dof_f_xmax; % cell and face on Neumann bnd
Param.qb= -q3(Param.dof_neu);  % vol. flux at Neumann bnd
g=[h1(0);h3(100)];% head a Dirichlet bnd
% Solve BLVP
[G,D,I] = build_ops(Grid);
L = -D*KdA*G;
fs =0;
[B,N,fn] = build_bnd(Param,Grid,I);
ha = solve_lbvp(L,fs+fn,B,g,N);
qa = comp_flux(D,KdA,G,ha,fs,Grid,Param);
%% Harmonic
L = -D*KdH*G;
fs =0;
[B,N,fn] = build_bnd(Param,Grid,I);
hh = solve_lbvp(L,fs+fn,B,g,N);
qh = comp_flux(D,KdH,G,hh,fs,Grid,Param);
%% Plotting...
subplot(1,2,1)
hold on
plot(Grid.xc,ha,'r.')
plot(Grid.xc,hh,'g*')
plot(Grid.xc,[h1(Grid.xc(1:P1)),h2(Grid.xc(P1+1:P2)),h3(Grid.xc(P2+1:end))]);
xlim([0 100]);
hold off
%%
subplot(1,2,2)
hold on 
plot(Grid.xf,[q1(0)*ones(1,P1+1),q2(0)*ones(1,P2+1-(P1+1)),q3(0)*ones(1,Grid.Nx-(P2))]);
plot(Grid.xf,qa,'r.','lineWidth',2)
plot(Grid.xf,qh,'g*','lineWidth',1)
hold off
