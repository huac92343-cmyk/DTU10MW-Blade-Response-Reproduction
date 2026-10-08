function loads = centrifugal_loads(data, r, mass_per_length, Omega)

arguments
    data (1,1) struct
    r double {mustBeNonnegative}
    mass_per_length double {mustBeNonnegative}
    Omega (1,1) double {mustBeNonnegative}
end

r = r(:); m = mass_per_length(:);
if numel(r) ~= numel(m)
    error('DTU10MW:SizeMismatch', 'r and mass_per_length must have equal lengths.');
end
radial = m.*Omega.^2.*r;
loads.flap = -radial*sin(data.cone_rad);
loads.edge = zeros(size(radial));
loads.span = radial*cos(data.cone_rad);
end
