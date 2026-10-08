function loads = gravity_loads(data, mass_per_length, azimuth_rad)

arguments
    data (1,1) struct
    mass_per_length double {mustBeNonnegative}
    azimuth_rad (1,1) double {mustBeFinite}
end

m = mass_per_length(:);
kin = coordinate_transform_blade(data, azimuth_rad, zeros(3,1));
gravity_wake = [0;0;-data.g];
gravity_blade = kin.R_rotor_to_blade * kin.R_wake_to_rotor * gravity_wake;
loads.flap = m*gravity_blade(1);
loads.edge = m*gravity_blade(2);
loads.span = m*gravity_blade(3);
loads.acceleration_blade = gravity_blade;
end
