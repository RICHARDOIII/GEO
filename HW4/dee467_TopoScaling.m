clc,clear all, close all;
L=[500 600 250 310 154 367 412 312 232 2032];
H=[50 60 60 45  47 39 36 86 119 119];
dh=[5 15 12 6 6 18 18 11 6 14];
Kh = [1e-6 5e-6 3e-7 4e-7 1e-7 1e-7 1e-7 1e-8 2e-8 2e-7];
Kv = [1e-6 5e-7 9e-9 8e-9 2e-9 5e-9 7e-9 3e-11 3e-11 3e-8];
ho=H;
a =(L./H).^2.*(Kv./Kh)
z=[2.8 5.65 39.75 24.5 36.75 6.367 3.55 77.9 114.5 7.25];
D=H-z;
disp('D/H')
DH=D./H;
disp(DH)
n=log(DH)./log(a);
loglog(DH,a.^n);
xlim([0 1])

