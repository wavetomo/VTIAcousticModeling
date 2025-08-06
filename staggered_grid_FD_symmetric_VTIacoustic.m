clc, clear;
% Description:
% Version: 1.0
% Autor: WaveTomo
% Date: 2022-12-31
% LastEditors: ZhangPingMin
% LastEditTime: 2025-08-01
% solving the 2D symmetric VTI acoustic with staggered finite difference

%% Marmousi Model
nz = 351;
nx = 1301;
nt = 2001;
dx = 10;
dz = 10;
dt = 1e-3;

vp = read_matrix('./MODEL_P-WAVE_VELOCITY_351x1301_10m.vel', nz, nx);
rho = read_matrix('./MODEL_DENSITY_351x1301_10m.vel', nz, nx);
delta = read_matrix('./MODEL_DELTA_351x1301_10m.vel', nz, nx);
epsilon = read_matrix('./MODEL_EPSLION_351x1301_10m.vel', nz, nx);

%% Hess Model
% nz = 366;
% nx = 882;
% nt = 3001;
% dx = 25;
% dz = 25;
% dt = 1.5e-3;
%
% vp = read_matrix('./vp_hess_366_882_25m.bin', nz, nx);
% rho = read_matrix('./rho_hess_366_882_25m.bin', nz, nx);
% delta = read_matrix('./delta_hess_366_882_25m.bin', nz, nx);
% epsilon = read_matrix('./epsilon_hess_366_882_25m.bin', nz, nx);

%% constant Model
% nz = 351;
% nx = 1301;
% nt = 1601;
% dx = 10;
% dz = 10;
% dt = 1e-3;
%
% vp = ones(nz, nx)*2000;
% rho = ones(nz, nx)*1500;
% delta = ones(nz, nx)*0.1;
% epsilon = ones(nz, nx)*0.2;

%%
x = (0:nx - 1) * dx;
z = (0:nz - 1) * dz;
t = (0:nt - 1) * dt;

nlayer = 30;
Nz = nz + 2 * nlayer;
Nx = nx + 2 * nlayer;

fpeak = 10;
src_z = 5;
% src_z = nz - 5;
src_x = round(nx/2);

rec_z = 1;

dsnap = 200;
nsnap = length(1:dsnap:nt);
fileName = 'VTIacoustic2d02.mat';

%%
src_z = src_z + nlayer;
src_x = src_x + nlayer;
rec_z = rec_z + nlayer;
src = Ricker(t, fpeak);

%%
vp = modpad2d(vp, nlayer, Nz, Nx);
rho = modpad2d(rho, nlayer, Nz, Nx);
delta = modpad2d(delta, nlayer, Nz, Nx);
epsilon = modpad2d(epsilon, nlayer, Nz, Nx);
[D, Beta] = meal2d(vp, nlayer, dx, dz, fpeak, Nx, Nz);

D = D ./ Beta;
rho = rho .* Beta;
vp = vp ./ Beta;

coef1 = (1 - (D * dt) / 2) ./ (1 + D * dt / 2);
coef2 = dt ./ (1 + D * dt / 2);

%%
invrhox = zeros(Nz, Nx);
invrhoz = zeros(Nz, Nx);
Kepsilon = rho .* (2 * epsilon) .* vp.^2;
Kdelta1 = 0.5 * (rho .* (1.0 + sqrt(1+2*delta)) .* vp.^2);
Kdelta2 = 0.5 * (rho .* (1.0 - sqrt(1+2*delta)) .* vp.^2);
for ix = 1:Nx
    for iz = 1:Nz
        if ix < Nx
            invrhox(iz, ix) = 2.0 / (rho(iz, ix) + rho(iz, ix+1));
        else
            invrhox(iz, ix) = 1.0 / rho(iz, ix);
        end
        if iz < Nz
            invrhoz(iz, ix) = 2.0 / (rho(iz, ix) + rho(iz+1, ix));
        else
            invrhoz(iz, ix) = 1.0 / rho(iz, ix);
        end
    end
end

%% 10th-order finite-differnce
N = 5;
c1 = 1.236425;
c2 = -0.1081130;
c3 = 0.02339911;
c4 = -0.5061550e-2;
c5 = 0.7054313e-3;
invdx = 1 / dx;
invdz = 1 / dz;

%%
vx = zeros(Nz, Nx);
vz = zeros(Nz, Nx);
p = zeros(Nz, Nx);
q = zeros(Nz, Nx);
r = zeros(Nz, Nx);
sigma_xx = zeros(Nz, Nx);
sigma_zz = zeros(Nz, Nx);

snaps = zeros(nz, nx, nsnap);
seis = zeros(nt, nx);

%%
k = 0;
izStart = N + 1;
izEnd = Nz - N;
ixStart = N + 1;
ixEnd = Nx - N;

for it = 1:nt
    for ix = ixStart:ixEnd
        for iz = izStart:izEnd
            dsigma_xxdx = c1 * (sigma_xx(iz, ix+1) - sigma_xx(iz, ix)) + ...
                c2 * (sigma_xx(iz, ix+2) - sigma_xx(iz, ix-1)) + ...
                c3 * (sigma_xx(iz, ix+3) - sigma_xx(iz, ix-2)) + ...
                c4 * (sigma_xx(iz, ix+4) - sigma_xx(iz, ix-3)) + ...
                c5 * (sigma_xx(iz, ix+5) - sigma_xx(iz, ix-4));
            dsigma_zzdz = c1 * (sigma_zz(iz+1, ix) - sigma_zz(iz, ix)) + ...
                c2 * (sigma_zz(iz+2, ix) - sigma_zz(iz-1, ix)) + ...
                c3 * (sigma_zz(iz+3, ix) - sigma_zz(iz-2, ix)) + ...
                c4 * (sigma_zz(iz+4, ix) - sigma_zz(iz-3, ix)) + ...
                c5 * (sigma_zz(iz+5, ix) - sigma_zz(iz-4, ix));

            dsigma_xxdx = dsigma_xxdx * invdx;
            dsigma_zzdz = dsigma_zzdz * invdz;

            vx(iz, ix) = coef1(iz, ix) * vx(iz, ix) + coef2(iz, ix) * ...
                invrhox(iz, ix) * dsigma_xxdx;
            vz(iz, ix) = coef1(iz, ix) * vz(iz, ix) + coef2(iz, ix) * ...
                invrhoz(iz, ix) * dsigma_zzdz;
        end
    end
    for ix = ixStart:ixEnd
        for iz = izStart:izEnd
            dvxdx = c1 * (vx(iz, ix) - vx(iz, ix-1)) + ...
                c2 * (vx(iz, ix+1) - vx(iz, ix-2)) + ...
                c3 * (vx(iz, ix+2) - vx(iz, ix-3)) + ...
                c4 * (vx(iz, ix+3) - vx(iz, ix-4)) + ...
                c5 * (vx(iz, ix+4) - vx(iz, ix-5));
            dvzdz = c1 * (vz(iz, ix) - vz(iz-1, ix)) + ...
                c2 * (vz(iz+1, ix) - vz(iz-2, ix)) + ...
                c3 * (vz(iz+2, ix) - vz(iz-3, ix)) + ...
                c4 * (vz(iz+3, ix) - vz(iz-4, ix)) + ...
                c5 * (vz(iz+4, ix) - vz(iz-5, ix));

            dvxdx = dvxdx * invdx;
            dvzdz = dvzdz * invdz;

            p(iz, ix) = coef1(iz, ix) * p(iz, ix) + coef2(iz, ix) * ...
                Kdelta1(iz, ix) * (dvxdx + dvzdz);
            q(iz, ix) = coef1(iz, ix) * q(iz, ix) + coef2(iz, ix) * ...
                Kdelta2(iz, ix) * (dvxdx - dvzdz);
            r(iz, ix) = coef1(iz, ix) * r(iz, ix) + coef2(iz, ix) * ...
                Kepsilon(iz, ix) * dvxdx;
        end
    end
    % load source
    p(src_z, src_x) = p(src_z, src_x) + 0.5 * src(it) * invdx * invdz * dt;
    % convert wavfield
    sigma_xx = p + q + r;
    sigma_zz = p - q;
    % reciever
    seis(it, :) = 0.5 * (sigma_xx(rec_z, nlayer+1:nlayer+nx) + sigma_zz(rec_z, nlayer+1:nlayer+nx));
    % snapshot
    if mod(it-1, dsnap) == 0 && it > 1
        k = k + 1;
        snap = 0.5 * (boundary2dcut(sigma_xx, nz, nx, nlayer) + ...
            boundary2dcut(sigma_zz, nz, nx, nlayer));
        snaps(:, :, k) = snap;

        clip = 0.5 * max(abs(snap(:)));
        figure(1);
        imagesc(x/1000, z/1000, snap), axis image;
        colormap(gca, gray), colorbar, clim([-1, 1]*clip);
        title(['Symmetric: P Snapshot at ', num2str(dt*(it - 1)*1000), ' ms']);
        xlabel('Distance / km'), ylabel('Depth / km');
    end
end

%%
% save(fileName, 'snap', 'seis');