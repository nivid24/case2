% Parameters
Lx = 1;  % Length of the plate in x-direction
Ly = 1;  % Length of the plate in y-direction
Lz = 1;  % Length of the plate in z-direction
Nx = 20; % Number of grid points in x-direction
Ny = 20; % Number of grid points in y-direction
Nz = 20; % Number of grid points in z-direction
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dz = Lz / (Nz - 1);
alpha = 0.1; % Thermal diffusivity
dt = 0.001;  % Time step size
timesteps = 1000;  % Number of time steps

% Initial condition
T = zeros(Nx, Ny, Nz);
T(:, :, :) = 100;  % Initial temperature

% Adding a heat flux at the center
center_x = floor(Nx / 2);
center_y = floor(Ny / 2);
center_z = floor(Nz / 2);
heat_flux = 100; % Magnitude of the heat flux

% Boundary conditions
T(:, :, 1) = 0;     % Bottom boundary (Dirichlet)
T(:, :, end) = 0;   % Top boundary (Dirichlet)
T(1, :, :) = 50;    % Left boundary (Dirichlet)
T(end, :, :) = 50;  % Right boundary (Dirichlet)
T(:, 1, :) = 50;    % Front boundary (Dirichlet)
T(:, end, :) = 50;  % Back boundary (Dirichlet)

% Finite difference method
for t = 1:timesteps
    T_old = T;
    % Adding heat flux at the center
    T(center_x, center_y, center_z) = T_old(center_x, center_y, center_z) + heat_flux * dt / (alpha * dx^2);
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 2:Nz-1
                T(i, j, k) = T_old(i, j, k) + alpha * dt * ((T_old(i+1, j, k) - 2*T_old(i, j, k) + T_old(i-1, j, k))/dx^2 + ...
                    (T_old(i, j+1, k) - 2*T_old(i, j, k) + T_old(i, j-1, k))/dy^2 + ...
                    (T_old(i, j, k+1) - 2*T_old(i, j, k) + T_old(i, j, k-1))/dz^2);
            end
        end
    end
end

% Plotting temperature colormap
[X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
figure;
slice(X, Y, Z, T, Nx/2, Ny/2, Nz/2);
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar;
title('Temperature Colormap with Heat Flux at the Center');