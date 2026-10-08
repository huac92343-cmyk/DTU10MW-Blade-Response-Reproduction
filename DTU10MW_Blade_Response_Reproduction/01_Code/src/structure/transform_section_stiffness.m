function stiffness = transform_section_stiffness(data, r, pitch_rad, ...
    elastic_twist_rad)
% EI_flap_Nm2 and EI_edge_Nm2 are the supplied principal bending stiffnesses.
% GJ is also supplied input data. Loads do not determine these base properties.
% Angles passed to this function are in radians; returned stiffnesses are in N*m^2.
% The section angle is geometric twist minus structural pitch plus pitch and elastic twist.

arguments
    data (1,1) struct
    r double
    pitch_rad (1,1) double {mustBeFinite}
    elastic_twist_rad double
end

r = r(:);
if isscalar(elastic_twist_rad)
    elastic_twist_rad = repmat(elastic_twist_rad,size(r));
else
    elastic_twist_rad = elastic_twist_rad(:);
end
if numel(elastic_twist_rad) ~= numel(r)
    error('DTU10MW:TwistSize', ...
        'elastic_twist_rad must be scalar or have one value per requested station.');
end

s = data.structure;
EI_x = interp1(s.r_m,s.EI_flap_Nm2,r,'pchip');
EI_y = interp1(s.r_m,s.EI_edge_Nm2,r,'pchip');
GJ = interp1(s.r_m,s.GJ_Nm2,r,'pchip');
structural_pitch = deg2rad(interp1(s.r_m,s.structural_pitch_deg,r,'pchip'));
geometric_twist = deg2rad(interp1(data.geometry.r_m, ...
    data.geometry.aerodynamic_twist_deg,r,'pchip'));

baseline_axis = geometric_twist-structural_pitch;
gamma_increment = pitch_rad+elastic_twist_rad;
total_axis = baseline_axis+gamma_increment;
c = cos(total_axis);
sine = sin(total_axis);

stiffness.r_m = r;
stiffness.EI_principal_x_Nm2 = EI_x;
stiffness.EI_principal_y_Nm2 = EI_y;
stiffness.GJ_Nm2 = GJ;
stiffness.geometric_twist_rad = geometric_twist;
stiffness.structural_pitch_rad = structural_pitch;
stiffness.baseline_principal_axis_rad = baseline_axis;
stiffness.pitch_rad = repmat(pitch_rad,size(r));
stiffness.elastic_twist_rad = elastic_twist_rad;
stiffness.gamma_increment_rad = gamma_increment;
stiffness.total_axis_angle_rad = total_axis;
stiffness.EI_flap_effective_Nm2 = EI_x.*c.^2+EI_y.*sine.^2;
stiffness.EI_edge_effective_Nm2 = EI_x.*sine.^2+EI_y.*c.^2;
stiffness.EI_flap_edge_coupling_Nm2 = (EI_x-EI_y).*sine.*c;
end
