function twist = compute_elastic_twist(data, bem)
% r is the packaged radial coordinate [m]; GJ is torsional stiffness [N*m^2].
% Direct torsion contains the CM moment and aerodynamic force-arm moment.
% Total torsion also contains the shear-center term multiplied by the configured factor.
% Internal moments are integrated from tip to root; twist is integrated from root to tip.
% psi_structure_rad and psi_aero_rad contain direct elastic twist for feedback.
% elastic_twist_total_rad contains total elastic twist, excluding geometric twist and pitch.

arguments
    data (1,1) struct
    bem (1,1) struct
end

required_structure_fields = {'GJ_Nm2','GA_flap_N','GA_edge_N', ...
    'shear_center_x_m','shear_center_y_m'};
if ~all(ismember(required_structure_fields, ...
        data.structure.Properties.VariableNames))
    error('DTU10MW:TwistCouplingData', ...
        'The structural GJ, GA and shear-center fields are required.');
end
if ~isfield(data,'twist_coupling_factor') || ...
        ~isscalar(data.twist_coupling_factor) || ...
        ~isfinite(data.twist_coupling_factor)
    error('DTU10MW:TwistCouplingFactor', ...
        'A finite scalar twist-coupling factor is required.');
end

r_aero = bem.r(:);
ea = interp1(data.structure.r_m, ...
    data.structure.elastic_axis_aft_le_chord,r_aero,'pchip');
if any(ea < 0 | ea > 1)
    error('DTU10MW:AxisLocation', ...
        'Elastic-axis locations must remain inside [0,1] chord.');
end
ac = interp1(data.structure.r_m, ...
    data.structure.aero_reference_axis_aft_le_chord,r_aero,'pchip');
arm_m = (ac-ea).*bem.chord;
distributed_torque_cm = -bem.pitching_moment_ref;
distributed_torque_force_arm = arm_m.*bem.normal_load;
r = double(data.structure.r_m(:));
distributed_torque_cm_structure = interp1(r_aero,distributed_torque_cm,r,'pchip');
distributed_torque_force_arm_structure = interp1( ...
    r_aero,distributed_torque_force_arm,r,'pchip');
distributed_torque_direct = distributed_torque_cm_structure+ ...
    distributed_torque_force_arm_structure;
internal_torque_direct = integrate_from_tip(r,distributed_torque_direct);
internal_torque_cm = integrate_from_tip(r,distributed_torque_cm_structure);
internal_torque_force_arm = integrate_from_tip( ...
    r,distributed_torque_force_arm_structure);

normal_load_structure = interp1(r_aero,bem.normal_load,r,'pchip');
tangential_load_structure = interp1(r_aero,bem.tangential_load,r,'pchip');
internal_shear_flap = integrate_from_tip(r,normal_load_structure);
internal_shear_edge = integrate_from_tip(r,tangential_load_structure);
shear_center_x = double(data.structure.shear_center_x_m(:));
shear_center_y = double(data.structure.shear_center_y_m(:));
internal_torque_coupling_raw = -shear_center_x.*internal_shear_flap+ ...
    shear_center_y.*internal_shear_edge;
coupling_factor = double(data.twist_coupling_factor);
internal_torque_coupling_corrected = coupling_factor.* ...
    internal_torque_coupling_raw;
internal_torque_total = internal_torque_direct+ ...
    internal_torque_coupling_corrected;

GJ = double(data.structure.GJ_Nm2(:));
kappa_direct = internal_torque_direct./GJ;
kappa_coupling_raw = internal_torque_coupling_raw./GJ;
kappa_coupling_corrected = coupling_factor.*kappa_coupling_raw;
kappa_total = internal_torque_total./GJ;
psi_direct = cumtrapz(r,kappa_direct);
psi_coupling_raw = cumtrapz(r,kappa_coupling_raw);
psi_coupling_corrected = cumtrapz(r,kappa_coupling_corrected);
psi_total = cumtrapz(r,kappa_total);
psi_cm = cumtrapz(r,internal_torque_cm./GJ);
psi_force_arm = cumtrapz(r,internal_torque_force_arm./GJ);

twist.r_structure = r;
twist.normal_load_structure_Npm = normal_load_structure;
twist.tangential_load_structure_Npm = tangential_load_structure;
twist.internal_shear_flap_N = internal_shear_flap;
twist.internal_shear_edge_N = internal_shear_edge;
twist.shear_center_x_m = shear_center_x;
twist.shear_center_y_m = shear_center_y;
twist.distributed_torque_Nm_m = distributed_torque_direct;
twist.distributed_torque_cm_Nm_m = distributed_torque_cm_structure;
twist.distributed_torque_force_arm_Nm_m = ...
    distributed_torque_force_arm_structure;
twist.internal_torque_Nm = internal_torque_direct;
twist.internal_torque_cm_Nm = internal_torque_cm;
twist.internal_torque_force_arm_Nm = internal_torque_force_arm;
twist.internal_torque_direct_Nm = internal_torque_direct;
twist.internal_torque_coupling_Nm = internal_torque_coupling_raw;
twist.internal_torque_coupling_raw_Nm = internal_torque_coupling_raw;
twist.internal_torque_coupling_corrected_Nm = ...
    internal_torque_coupling_corrected;
twist.internal_torque_total_Nm = internal_torque_total;
twist.GJ_Nm2 = GJ;
twist.twist_coupling_factor = coupling_factor;
twist.kappa_twist_direct_radpm = kappa_direct;
twist.kappa_twist_coupling_radpm = kappa_coupling_raw;
twist.kappa_twist_coupling_raw_radpm = kappa_coupling_raw;
twist.kappa_twist_coupling_corrected_radpm = ...
    kappa_coupling_corrected;
twist.kappa_twist_total_radpm = kappa_total;
twist.elastic_twist_direct_rad = psi_direct;
twist.elastic_twist_coupling_rad = psi_coupling_raw;
twist.elastic_twist_coupling_raw_rad = psi_coupling_raw;
twist.elastic_twist_coupling_corrected_rad = psi_coupling_corrected;
twist.elastic_twist_total_rad = psi_total;
twist.psi_structure_rad = psi_direct;
twist.psi_cm_structure_rad = psi_cm;
twist.psi_force_arm_structure_rad = psi_force_arm;
twist.psi_aero_rad = interp1(r,psi_direct,r_aero,'pchip');
twist.psi_total_aero_rad = interp1(r,psi_total,r_aero,'pchip');
twist.elastic_axis_aft_le_chord = ea;
twist.aero_reference_axis_aft_le_chord = ac;
twist.force_arm_m = arm_m;
twist.tip_twist_rad = psi_direct(end);
twist.tip_twist_deg = rad2deg(psi_direct(end));
twist.tip_twist_direct_rad = psi_direct(end);
twist.tip_twist_coupling_rad = psi_coupling_raw(end);
twist.tip_twist_coupling_raw_rad = psi_coupling_raw(end);
twist.tip_twist_coupling_corrected_rad = psi_coupling_corrected(end);
twist.tip_twist_total_rad = psi_total(end);
twist.tip_twist_direct_deg = rad2deg(psi_direct(end));
twist.tip_twist_coupling_deg = rad2deg(psi_coupling_raw(end));
twist.tip_twist_coupling_raw_deg = rad2deg(psi_coupling_raw(end));
twist.tip_twist_coupling_corrected_deg = ...
    rad2deg(psi_coupling_corrected(end));
twist.tip_twist_total_deg = rad2deg(psi_total(end));
twist.tip_twist_cm_rad = psi_cm(end);
twist.tip_twist_force_arm_rad = psi_force_arm(end);
twist.tip_twist_cm_deg = rad2deg(psi_cm(end));
twist.tip_twist_force_arm_deg = rad2deg(psi_force_arm(end));
end

function internal = integrate_from_tip(r,distributed)
internal = zeros(size(r));
for j = numel(r)-1:-1:1
    internal(j) = internal(j+1)+0.5*(distributed(j)+distributed(j+1))* ...
        (r(j+1)-r(j));
end
end
