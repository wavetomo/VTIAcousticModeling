clc, clear;
% Description:
% Version: 1.0
% Autor: WaveTomo
% Date: 2022-12-31
% LastEditors: ZhangPingMin
% LastEditTime: 2025-08-01
% solving the 2D non symmetric VTI acoustic with staggered finite difference

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
fileName = 'VTIacoustic2d01.mat';

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
c11 = rho .* (1 + 2 * epsilon) .* vp.^2;
c13 = rho .* sqrt(1+2*delta) .* vp.^2;
c33 = rho .* vp.^2;
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

%% 差分系数
% 优化差分系数
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
sigma_xx = zeros(Nz, Nx);
sigma_zz = zeros(Nz, Nx);


Vx = zeros(Nz, Nx);
Vz = zeros(Nz, Nx);
tau_xx = zeros(Nz, Nx);
tau_zz = zeros(Nz, Nx);

snaps = zeros(nz, nx, nsnap);
seis = zeros(nt, nx);

%%
k = 0;
izStart = N + 1;
izEnd = Nz - N;
ixStart = N + 1;
ixEnd = Nx - N;

for it = 1:nt
    tau_xx = c11 .* sigma_xx + c13 .* sigma_zz;
    tau_zz = c13 .* sigma_xx + c33 .* sigma_zz;
    for ix = ixStart:ixEnd
        for iz = izStart:izEnd
            dtau_xxdx = c1 * (tau_xx(iz, ix+1) - tau_xx(iz, ix)) + ...
                c2 * (tau_xx(iz, ix+2) - tau_xx(iz, ix-1)) + ...
                c3 * (tau_xx(iz, ix+3) - tau_xx(iz, ix-2)) + ...
                c4 * (tau_xx(iz, ix+4) - tau_xx(iz, ix-3)) + ...
                c5 * (tau_xx(iz, ix+5) - tau_xx(iz, ix-4));
            dtau_zzdz = c1 * (tau_zz(iz+1, ix) - tau_zz(iz, ix)) + ...
                c2 * (tau_zz(iz+2, ix) - tau_zz(iz-1, ix)) + ...
                c3 * (tau_zz(iz+3, ix) - tau_zz(iz-2, ix)) + ...
                c4 * (tau_zz(iz+4, ix) - tau_zz(iz-3, ix)) + ...
                c5 * (tau_zz(iz+5, ix) - tau_zz(iz-4, ix));

            dtau_xxdx = dtau_xxdx * invdx;
            dtau_zzdz = dtau_zzdz * invdz;

            vx(iz, ix) = coef1(iz, ix) * vx(iz, ix) + coef2(iz, ix) * dtau_xxdx;
            vz(iz, ix) = coef1(iz, ix) * vz(iz, ix) + coef2(iz, ix) * dtau_zzdz;
        end
    end
    Vx = invrhox .* vx;
    Vz = invrhoz .* vz;
    for ix = ixStart:ixEnd
        for iz = izStart:izEnd
            dVxdx = c1 * (Vx(iz, ix) - Vx(iz, ix-1)) + ...
                c2 * (Vx(iz, ix+1) - Vx(iz, ix-2)) + ...
                c3 * (Vx(iz, ix+2) - Vx(iz, ix-3)) + ...
                c4 * (Vx(iz, ix+3) - Vx(iz, ix-4)) + ...
                c5 * (Vx(iz, ix+4) - Vx(iz, ix-5));
            dVzdz = c1 * (Vz(iz, ix) - Vz(iz-1, ix)) + ...
                c2 * (Vz(iz+1, ix) - Vz(iz-2, ix)) + ...
                c3 * (Vz(iz+2, ix) - Vz(iz-3, ix)) + ...
                c4 * (Vz(iz+3, ix) - Vz(iz-4, ix)) + ...
                c5 * (Vz(iz+4, ix) - Vz(iz-5, ix));

            dVxdx = dVxdx * invdx;
            dVzdz = dVzdz * invdz;

            sigma_xx(iz, ix) = coef1(iz, ix) * sigma_xx(iz, ix) + coef2(iz, ix) * dVxdx;
            sigma_zz(iz, ix) = coef1(iz, ix) * sigma_zz(iz, ix) + coef2(iz, ix) * dVzdz;
        end
    end
    % load source
    sigma_xx(src_z, src_x) = sigma_xx(src_z, src_x) + 0.5 * src(it) * invdx * invdz * dt;
    sigma_zz(src_z, src_x) = sigma_zz(src_z, src_x) + 0.5 * src(it) * invdx * invdz * dt;
    % reciever
    seis(it, :) = 0.5 * (sigma_xx(rec_z, nlayer+1:nlayer+nx) + sigma_zz(rec_z, nlayer+1:nlayer+nx));
    % snapshot
    if mod(it-1, dsnap) == 0 && it > 1
        k = k + 1;
        snap = 0.5 * (boundary2dcut(sigma_xx, nz, nx, nlayer) + ...
            boundary2dcut(sigma_zz, nz, nx, nlayer));
        snaps(:, :, k) = snap;

        clip = 0.5 * max(abs(snap(:)));
        figure(2);
        imagesc(x/1000, z/1000, snap), axis image;
        colormap(gca, gray), colorbar, clim([-1, 1]*clip);
        title(['Non-symmetric: P Snapshot at ', num2str(dt*(it - 1)*1000), ' ms']);
        xlabel('Distance / km'), ylabel('Depth / km');
    end

end

%%
% save(fileName, 'snap', 'seis');