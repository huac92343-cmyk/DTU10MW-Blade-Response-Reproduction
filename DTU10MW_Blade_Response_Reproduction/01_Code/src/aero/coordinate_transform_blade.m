function kin = coordinate_transform_blade(data, azimuth_rad, wind_vector_wake)

arguments
    data (1,1) struct
    azimuth_rad (1,1) double {mustBeFinite}
    wind_vector_wake (3,1) double {mustBeFinite}
end

cy = cos(data.yaw_rad);   sy = sin(data.yaw_rad);
ct = cos(data.shaft_tilt_rad); st = sin(data.shaft_tilt_rad);
R_yaw = [cy sy 0; -sy cy 0; 0 0 1];
R_tilt = [ct 0 -st; 0 1 0; st 0 ct];
wind_rotor = R_tilt * R_yaw * wind_vector_wake;

e_shaft = [1;0;0];
e_radial = [0; sin(azimuth_rad); cos(azimuth_rad)];
e_edge = [0; cos(azimuth_rad); -sin(azimuth_rad)];
cb = cos(data.cone_rad); sb = sin(data.cone_rad);
e_blade = sb*e_shaft + cb*e_radial;
e_flap = cb*e_shaft - sb*e_radial;

R_rotor_to_blade = [e_flap e_edge e_blade].';
orthogonality_error = norm(R_rotor_to_blade*R_rotor_to_blade.'-eye(3),'fro');
if orthogonality_error > 1e-12
    error('DTU10MW:CoordinateTransform', ...
        'Blade basis is not orthonormal (error %.3g).', orthogonality_error);
end

kin.R_wake_to_rotor = R_tilt * R_yaw;
kin.R_rotor_to_blade = R_rotor_to_blade;
kin.wind_rotor = wind_rotor;
kin.wind_blade = R_rotor_to_blade * wind_rotor;
kin.e_shaft = e_shaft;
kin.e_flap = e_flap;
kin.e_edge = e_edge;
kin.e_blade = e_blade;
kin.orthogonality_error = orthogonality_error;
end
