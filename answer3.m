% mohammad khalaji 40105504
% hw1---answer 3 a
clc, clear all

f = 2e9; % Frequency
c = 3e8; % Speed of light 
lambda_0 = c / f; % Wavelength 
d = 2 * lambda_0; % the slab Thickness
epsilon_r = -1; % in DNG
mu_r = -1; % in DNG
n_0 = 1; % Refractive index of air
n_1 = -sqrt(epsilon_r * mu_r); % Refractive index of DNG slab ( this negative )
eta_0 = 120 * pi;
eta_1 = 120 * pi * sqrt( mu_r/epsilon_r );

% angle of the incident plane wave
theta_in = 60 * pi / 180; 

% angle of the refraction olane wave
theta_out = asin( n_0 / n_1 * sin(theta_in) ); % snell's law
 



% generation meshgrid
x = linspace(-2 * lambda_0, 4 * lambda_0, 1000);
z = linspace(-2 * lambda_0, 4 * lambda_0, 1000);
[X, Z] = meshgrid(x, z);

% the wave vectors
k_0 = 2 * pi / lambda_0; % Wave number in free space
k_1 = 2 * pi / lambda_0; % wave number in DNG
kx_0 = k_0 * sin(theta_in); % Incident wave x-component
kz_0 = k_0 * cos(theta_in); % Incident wave z-component
kx_1 = k_0 * n_1 * abs(sin(theta_out)); % Transmitted wave x-component in slab
kz_1 = k_0 * n_1 * abs(cos(theta_out)); % Transmitted wave z-component in slab

% boundary condition 
syms E_10 E_20  T R
E_i0 = 100; % magnetiude of the incident plane wave in v/m
eq1 = E_10 * exp(-i* kz_1 * d) + E_20 * exp(i * kz_1 *d) == E_i0 * T * exp(-i * kz_0 * d);
eq2 = - cos(theta_out) * E_10 / eta_1  * exp(-i* kz_1 * d) + cos(theta_out) * E_20 / eta_1 * exp(i * kz_1 *d) == - cos(theta_in) * E_i0 / eta_0 * T * exp(-i * kz_0 * d);
eq3 = - cos(theta_out) * E_10 / eta_1 + cos(theta_out) * E_20 / eta_1 == cos(theta_in) * E_i0 * R / eta_0 - cos(theta_in) * E_i0 / eta_0 ;
eq4 = E_10 + E_20  == E_i0 *R + E_i0;
[E_10 E_20 T R] = solve([eq1 eq2 eq3 eq4], [E_10 E_20 T R] );
E_10 = double(subs(E_10));
E_20 = double(subs(E_20));
T = double(subs(T));
R= double(subs(R));

% fields befor slab
E_i = E_i0 * exp(-1i * (kx_0 * X + kz_0 * Z));
H_i =  E_i0 / eta_0 * exp(-1i * (kx_0 * X + kz_0 * Z)) ;
E_r = E_i0 * R * exp(-1i * (kx_0 * X - kz_0 * Z));
H_r = E_i0 * R / eta_0 * exp(-1i * (kx_0 * X - kz_0 * Z));

%fields in slab
E_1 = E_10 * exp(-i * (-kx_1 * X + kz_1*Z ));
H_1 = E_10 /eta_1 * exp(i * (-kx_1 * X + kz_1*Z )) ;
E_2 = E_20 * exp(-i * (kx_1 * X - kz_1*Z ));
H_2 = E_20 /eta_1 * exp(i * (kx_1 * X - kz_1*Z )) ;

% fields after slab (Transmitted fields)_ 
E_t = E_i0 * T * exp(-i * (kx_0 * X + kz_0 * Z));
H_t = E_i0 / eta_0 * exp(-i * (kx_1 * X + kz_1 * Z));

% Plot the electric fields
figure('NAME','ELECTRIC FIELDS');
subplot(3,1,1);
imagesc(x, z, real(E_i+E_r));
colorbar;
title('Electric Field (befor slab) ');
% axis tight;

% Plot the electric fields
subplot(3, 1, 2);
imagesc(x, z, real(E_1+E_2));
colorbar;
colormap("jet")
title('Electric Field (IN the DNG slab) ');
% axis tight;

% Plot the electric fields
subplot(3, 1, 3);
imagesc(x, z, real(E_t));
colorbar;
title('Electric Field (after slab) ');
axis tight;

% Plot the magnetic fields
figure('NAME','MAGNETIC FIELDS');
subplot(3,1,1);
imagesc(x, z, real(H_i+H_r));
colorbar;
title('Magnetic Field (befor slab) ');
% axis tight;

% Plot the magnetic fields
subplot(3, 1, 2);
imagesc(x, z, real(H_1+H_2));
colorbar;
colormap("jet")
title('Magnetic Field (IN the DNG slab) ');
% axis tight;

% Plot the Magnetic fields
subplot(3, 1, 3);
imagesc(x, z, real(H_t));
colorbar;
title('Magnetic Field (after slab) ');
% axis tight;

