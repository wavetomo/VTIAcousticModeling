function [D, Beta] = meal2d(vp, nlayer, dx, dz, fm, Nx, Nz)
% Description: Modified 2D effective absorbing boundary condition
% Version: 1.0
% INPUT
% vp: velocity of P wave
% nlayer: number of absorbing boundary layers
% dx: samping interval along x dimention
% dz: samping interval along z dimention
% fm: peak
% NX: number of sampling points along z dimention
% NZ: number of sampling points along z dimention
% OUTPUT
% D: attenuation damping coefficient
% Beta: deceleration damping coefficient
% reference
% Yao, G., Da Silva, N. V., & Wu, D. (2018). An effective absorbing layer
% for the boundary condition in acoustic seismic wave simulation. Journal
% of Geophysics and Engineering, 15(2), 495-511.
% ----------------
% Autor: Zhang PingMin
% Date: 2022-12-31
% LastEditors: ZhangPingMin
% LastEditTime: 2023-01-04
% Copyright (c) 2023 WaveTomo. All rights reserved.

%%
alpha = 0.05 * nlayer + 1;
vp_max = max(max(vp));
Lx = nlayer * dx;
Lz = nlayer * dz;
R = 10^(-(log10(nlayer) - 1) / log10(2) - 3);
d0x = -log(R) * (3 * vp_max) / (2 * Lx);
d0z = -log(R) * (3 * vp_max) / (2 * Lz);
beta0x = vp_max / sqrt(3) / (5 * dx * fm);
beta0z = vp_max / sqrt(3) / (5 * dz * fm);

%%
ddz = zeros(Nz, 1);
ddx = zeros(1, Nx);
betaz = ones(Nz, 1);
betax = ones(1, Nx);
for i = 1:nlayer
    w = exp(log(2)*(((nlayer - i) / nlayer)^alpha)) - 1;
    ddz(i) = d0z * w;
    ddz(Nz-i+1) = ddz(i);
    ddx(i) = d0x * w;
    ddx(Nx-i+1) = ddx(i);
    betaz(i) = 1 + (beta0z - 1) * w;
    betaz(Nz-i+1) = betaz(i);
    betax(i) = 1 + (beta0x - 1) * w;
    betax(Nx-i+1) = betax(i);
end

%%
D = sqrt(ddx.^2+ddz.^2);
Beta = sqrt((betax - 1).^2+(betaz - 1).^2) + 1;
end