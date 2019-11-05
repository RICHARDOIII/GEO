clc,clear all, close all;
L=500;
H=50;dh=5;
Kh=1e-6;
Kv=Kh;
ho=H;
x = linspace(0,1,50);
[xd,zd]=meshgrid(x,x);
x=xd*L;
z=zd*H;
a =(L/H)^2*(Kv/Kh);
hd=@(xd,zd)((cos(2*pi*xd).*cosh(2*pi*zd./sqrt(a)))./cosh(2*pi/sqrt(a)));
Wd=@(xd,zd) sqrt(a)*((sin(2*pi*xd).*sinh(2*pi*zd./sqrt(a)))./cosh(2*pi/sqrt(a)));
%% Dimensionless what ever that means
hold on
h=@(a,b)hd(a,b)*dh+ho;
w=@(a,b)Wd(a,b)*Kh*dh*H/L;
contour(x,z,w(xd,zd),18);
hold on
contour(x,z,h(xd,zd),18);
