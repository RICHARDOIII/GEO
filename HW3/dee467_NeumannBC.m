% file name: eid_NeumannBC.m
% author: Pancho Villa
% date: 15 Apr 1915
clc, close all, clear all
%% Physical properties in cgs units
Length = 45;    % length of column [cm]
r = 2.5;         % radius of column [cm]
Q = 10;         % injection rate [cm^3/s] - note 1 mL = 1cm^3
phi = .4;       % porosity [1]
d = .023;         % bead diameter [cm]
rho = 1.0;         % water density [g/cm^3]
grav = 1;      % gravitational acceleration [gal = cm/s^2]
mu =  .01;       % water dyn. viscosity [poise = g/cm/s]
% Hydraulic conductivity from Kozeny-Carman
K = (rho*grav/mu)*(phi^3/((1-phi)^2))*(d^2/180);         % [cm/s]

%% Volumetric flux at boundary

A = r^2*pi;         % cross sectional area [cm^2]
qb = Q/A;        % [cm/s]

%%  Analytic solution

xa =0:.01:Length;         % fine grid for plot analytic solution
ha = @(x)-Q*x/(pi*r^2*K)+Q*Length/(pi*r^2*K);    % analytic vol. flux
qa = @(x)Q./(pi.*r.^2) ;    % analytic head%% Finite volume solution
Grid.xmin = 0; 
Grid.xmax = Length; 
Grid.Nx = 20;
Grid = build_grid(Grid);
%% Boundary conditions
Param.dof_dir =Grid.dof_xmin; 
Param.dof_f_dir =Grid.dof_f_xmin; 
% cell and face on Dirichlet bnd
Param.dof_neu=Grid.dof_xmax;
Param.dof_f_neu=Grid.dof_f_xmax; % cell and face on Neumann bnd
Param.qb= -qa(Param.dof_neu);  % vol. flux at Neumann bnd
g=[ha(Grid.xc(Grid.dof_xmin));ha(Grid.xc(Grid.dof_xmax))];% head a Dirichlet bnd
% Solve BLVP
[G,D,I] = build_ops(Grid);
L = -D*G;
fs =0;
[B,N,fn] = build_bnd(Param,Grid,I);
h = solve_lbvp(L,fs+fn,B,g,N);
q = comp_flux(D,K,G,h,fs,Grid,Param);

%% Plotting...
subplot(1,3,1)
hold on
plot(Grid.xc,h,'.')
plot(Grid.xc,ha(Grid.xc))
xlim([0 Length]);
hold off
subplot(1,3,2)
hold on 
fplot(qa);
ylim([.3 .6]);
xlim([0 Length]);
plot(Grid.xf,q,'.','lineWidth',2)
hold off
subplot(1,3,3)
hold on
y=Q;
line([0,Length],[y,y])
plot(Grid.xf,A*q,'r.')
ylim([8 11]);
xlim([0 Length]);
