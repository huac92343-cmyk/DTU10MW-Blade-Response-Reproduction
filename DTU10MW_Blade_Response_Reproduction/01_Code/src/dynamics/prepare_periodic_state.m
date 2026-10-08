function periodic = prepare_periodic_state(cfg, data, wind_speed, Omega, pitch_rad)

arguments
    cfg (1,1) struct
    data (1,1) struct
    wind_speed (1,1) double {mustBeNonnegative}
    Omega (1,1) double {mustBePositive}
    pitch_rad (1,1) double {mustBeFinite}
end

n_az = cfg.time.azimuth_samples;
azimuth = (0:n_az-1).'*2*pi/n_az;
states = cell(n_az,1);
psi = [];
for k = 1:n_az
    states{k} = solve_aeroelastic_equilibrium(cfg,data,wind_speed,Omega, ...
        pitch_rad,azimuth(k),psi);
    psi = states{k}.twist.psi_aero_rad;
end

ndof = numel(states{1}.force.free);
force_samples = zeros(ndof,n_az);
K = zeros(ndof);
C_aero = zeros(ndof);
tip_static = zeros(n_az,2);
tip_twist = zeros(n_az,1);
tip_twist_cm = zeros(n_az,1);
tip_twist_force_arm = zeros(n_az,1);
twist_structure = zeros(height(data.structure),n_az);
thrust = zeros(n_az,1);
torque = zeros(n_az,1);
for k = 1:n_az
    state = states{k};
    force_samples(:,k) = state.force.free;
    K = K+full(state.model.K)/n_az;
    aero_k = assemble_aerodynamic_damping(data,state.model,state.bem);
    C_aero = C_aero+full(aero_k.C)/n_az;
    tip_static(k,:) = [state.tip_flap_static_m state.tip_edge_static_m];
    tip_twist(k) = state.twist.tip_twist_rad;
    tip_twist_cm(k) = state.twist.tip_twist_cm_rad;
    tip_twist_force_arm(k) = state.twist.tip_twist_force_arm_rad;
    twist_structure(:,k) = state.twist.psi_structure_rad;
    thrust(k) = trapz(state.bem.r,state.bem.normal_load);
    torque(k) = trapz(state.bem.r,state.bem.tangential_load.*state.bem.r);
end

model = states{1}.model;
model.K = sparse(0.5*(K+K.'));
model.K_full(model.mesh.free_dofs,model.mesh.free_dofs) = model.K;
modes = solve_blade_modes(model,max(cfg.structure.mode_count, ...
    max(cfg.damping.rayleigh_mode_indices)));
aero_mean = aero_k;
aero_mean.C = sparse(0.5*(C_aero+C_aero.'));
damping = build_damping_matrix(cfg,data,model,modes,aero_mean);

periodic.wind_speed_m_s = wind_speed;
periodic.Omega_rad_s = Omega;
periodic.pitch_rad = pitch_rad;
periodic.azimuth_rad = azimuth;
periodic.force_samples = force_samples;
periodic.mean_force = mean(force_samples,2);
periodic.model = model;
periodic.modes = modes;
periodic.damping = damping;
periodic.tip_static_samples_m = tip_static;
periodic.tip_twist_samples_rad = tip_twist;
periodic.tip_twist_cm_samples_rad = tip_twist_cm;
periodic.tip_twist_force_arm_samples_rad = tip_twist_force_arm;
periodic.twist_structure_samples_rad = twist_structure;
tip_radius_m = double(data.R);
periodic.geometric_tip_twist_deg = interp1( ...
    double(data.geometry.r_m), ...
    double(data.geometry.aerodynamic_twist_deg),tip_radius_m,'pchip');
periodic.structural_tip_pitch_deg = interp1( ...
    double(data.structure.r_m), ...
    double(data.structure.structural_pitch_deg),tip_radius_m,'pchip');
periodic.baseline_principal_axis_tip_deg = ...
    periodic.geometric_tip_twist_deg-periodic.structural_tip_pitch_deg;
periodic.blade_thrust_samples_N = thrust;
periodic.blade_torque_samples_Nm = torque;
periodic.reference_state = states{1};
periodic.coupling_iterations = cellfun(@(s) s.coupling_iterations,states);
end
