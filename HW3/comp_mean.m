function [Kd] = comp_mean(K,p,Grid)
% author: Daniel Espinoza
% date: 14-Oct-2016
% Description:
% Takes coefficient field, K, defined at the cell centers and computes the
% mean specified by the power, p and returns it in a sparse diagonal
% matrix, Kd.

% p = power of the generalized mean
% 1 (arithmetic mean)
% -1 (harmonic mean)
if (p == -1) || (p == 1)
    if (Grid.Nx == Grid.N) || (Grid.Ny == Grid.N) % 1D
        mean = zeros(Grid.Nx+1,1);
        mean(2:Grid.Nx) =((K(1:end-1).^p+(K(2:end).^p))/2).^1/p;
        Kd = spdiags(mean,0,Grid.Nx+1,Grid.Nfx);
    elseif (Grid.N > Grid.Nx) || (Grid.N > Grid.Ny) % 2D
        error('2D permeability is not implemented')
    else
        error('3D permeability is not implemented')
    end
else
    error('This power does not have significance.')
end