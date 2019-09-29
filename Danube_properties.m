%% Physical properties
cm2m = 1/100;        % cm to m conversion
yr2s = 365*24*60^2;  % yr to s conversion

Length = 85070;      % Distance between Danube and Tisza rivers [m]
Width = 5430;        % Width of segment considered [m]
Param.K =2e-2*cm2m;        % Hydraulic conductivity [m/s]
qp = 1.5*cm2m/yr2s;  % Average annual precipitation [m3/m2/s]
Param.hL = 90;             % Elevation of Danube river[m]
Param.hR = 80;             % Elevation of Tisza river [m]
b = 100;             % Aquifer thickness [m]

% Data
mm2m = 1.81e3; % mm on map to m in reality
% dist_topo is the location of the topography measurements
dist_topo = [ 0 14  17  19  21 23  25  28  32  33  35  37  42 43 44 47]*mm2m;
% topo is the elevation of the land surface above sea-level
topo = [90 95 100 115 115 110 115 120 115 110 105 100 95 90 85 80];
% diest_gw is the location of the groundwater table measurements
dist_gw = [ 0 15 18 21 25 28 34 36 39 43 44 47]*mm2m;
% gw is the elevation of the gorundwater table above sealevel.
gw = [90 95 100 105 110 110 105 100 95 90 85 80];
fs=Param.qp/(b*Param.K);
Grid.xmin=0;Grid.xmax=Length;Grid.Nx=15;
Grid=build_grid(Grid);
[ G, D, I ] = build_ops( Grid );
L=-D*G;
Param.dof_dir = Grid.dof_xmin; % identify cells on Dirichlet bnd
Param.dof_f_dir = Grid.dof_f_xmin; % identify faces on Dirichlet bnd
Param.dof_neu = Grid.dof_xmax; % identify cells on Neumann bnd
Param.dof_f_neu = Grid.dof_f_xmax;
Param.qb=1;
Param.qp=qp;
Param.g=[Param.hL;Param.hR];
[B,N,fn] = build_bnd(Param,Grid,I);
h = solve_lbvp(L,fs+fn,B,Param.g,N);
% to plot the topography
% plot(dist_topo,topo,'-')