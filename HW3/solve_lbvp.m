function [u] = solve_lbvp(L,f,B,g,N)
% author: Daniel Espinoza
% date: 28/09/2019
% Description
% Computes the solution $u$ to the linear differential problem given by
hp= B'*((B*B')\g);
ho = N*((N.'*L*N)\N.'*(f-L*hp));
u = ho + hp;
end

