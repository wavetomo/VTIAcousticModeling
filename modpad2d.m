function [mod] = modpad2d(mod, nlayer, NZ, NX)
% Description: pading 2D model parameter for absorbing boundary condition
% Version: 1.0
% INPUT
% mod: model parameter
% nlayer: number of absorbing boundary layers
% NX: number of sampling points along z dimention
% NZ: number of sampling points along z dimention
% OUTPUT
% mod: model parameter after pading boundary
% ----------------
% Autor: Zhang PingMin
% Date: 2022-12-31
% LastEditors: ZhangPingMin
% LastEditTime: 2023-01-04
% Copyright (c) 2023 WaveTomo. All rights reserved.
[nz, nx] = size(mod);
modp = zeros(NZ, NX);
modp(nlayer+1:nz+nlayer, nlayer+1:nx+nlayer) = mod;
modp(1:nlayer, :) = modp(nlayer+1, :) .* ones(nlayer, NX);
modp(nz+nlayer+1:NZ, :) = modp(nz+nlayer, :) .* ones(nlayer, NX);
modp(:, 1:nlayer) = modp(:, nlayer+1) .* ones(NZ, nlayer);
modp(:, nx+nlayer+1:NX) = modp(:, nx+nlayer) .* ones(NZ, nlayer);
mod = modp;

end