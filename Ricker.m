function src = Ricker(t, fpeak, t0, A)

if nargin < 2
    fpeak = 10.0;
end

if nargin < 3
    t0 = 1.5 / fpeak;
end
if nargin < 4
    A = 1;
end

t = t - t0;

src = A * (1.0 - 2.0 * (pi * fpeak * t).^2.0) .* exp(-1.0*(pi * fpeak * t).^2.0);

end