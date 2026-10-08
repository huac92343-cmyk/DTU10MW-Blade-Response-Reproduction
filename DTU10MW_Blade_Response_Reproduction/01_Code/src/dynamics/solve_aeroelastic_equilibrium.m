function state = solve_aeroelastic_equilibrium(cfg, data, wind_speed, ...
    Omega, pitch_rad, azimuth_rad, initial_twist_aero_rad)

arguments
    cfg (1,1) struct
    data (1,1) struct
    wind_speed (1,1) double {mustBeNonnegative}
    Omega (1,1) double {mustBeNonnegative}
    pitch_rad (1,1) double {mustBeFinite}
    azimuth_rad (1,1) double {mustBeFinite}
    initial_twist_aero_rad double = []
end

n_aero = height(data.aerodynamic_stations);
if isempty(initial_twist_aero_rad)
    psi = zeros(n_aero,1);
else
    psi = initial_twist_aero_rad(:);
    if numel(psi) ~= n_aero
        error('DTU10MW:TwistSize', 'Initial twist must match aerodynamic stations.');
    end
end
previous_normal = [];
converged = false;
angle_error = inf;
load_error = inf;
for iteration = 1:cfg.coupling.max_iterations
    bem = bem_solver_DTU10MW(cfg,data,wind_speed,Omega,pitch_rad, ...
        azimuth_rad,psi);
    twist = compute_elastic_twist(data,bem);
    target = twist.psi_aero_rad;
    angle_error = max(abs(target-psi));
    if isempty(previous_normal)
        load_error = inf;
    else
        load_error = norm(bem.normal_load-previous_normal)/ ...
            max(norm(previous_normal),1);
    end
    psi = psi+cfg.coupling.relaxation*(target-psi);
    if angle_error <= cfg.coupling.tolerance_rad && ...
            load_error <= cfg.coupling.load_relative_tolerance
        converged = true;
        break
    end
    previous_normal = bem.normal_load;
end
if ~converged && cfg.coupling.fail_on_nonconvergence
    error('DTU10MW:AeroelasticNonconvergence', ...
        ['Aeroelastic twist iteration failed at V=%g m/s, azimuth=%g deg: ' ...
        'angle %.3g rad, load %.3g.'],wind_speed,rad2deg(azimuth_rad), ...
        angle_error,load_error);
end

bem = bem_solver_DTU10MW(cfg,data,wind_speed,Omega,pitch_rad, ...
    azimuth_rad,psi);
twist = compute_elastic_twist(data,bem);
model = assemble_blade_matrices(cfg,data,Omega,pitch_rad, ...
    twist.psi_structure_rad);
r_structure = double(data.structure.r_m(:));
mass = double(data.structure.mass_per_length_kg_m(:));
gravity = gravity_loads(data,mass,azimuth_rad);
centrifugal = centrifugal_loads(data,r_structure,mass,Omega);
aero_flap = interp1(bem.r,bem.normal_load,r_structure,'pchip');
aero_edge = interp1(bem.r,bem.tangential_load,r_structure,'pchip');
q_flap = aero_flap+gravity.flap+centrifugal.flap;
q_edge = aero_edge+gravity.edge+centrifugal.edge;
force = assemble_distributed_force(model,r_structure,q_flap,q_edge);
q_static_free = model.K\force.free;
q_static_full = zeros(model.mesh.ndof,1);
q_static_full(model.mesh.free_dofs) = q_static_free;

state.wind_speed_m_s = wind_speed;
state.Omega_rad_s = Omega;
state.pitch_rad = pitch_rad;
state.azimuth_rad = azimuth_rad;
state.bem = bem;
state.twist = twist;
state.model = model;
state.gravity = gravity;
state.centrifugal = centrifugal;
state.aero_flap_N_m = aero_flap;
state.aero_edge_N_m = aero_edge;
state.total_flap_N_m = q_flap;
state.total_edge_N_m = q_edge;
state.force = force;
state.q_static_full = q_static_full;
state.tip_flap_static_m = q_static_full(model.mesh.tip_dofs(1));
state.tip_edge_static_m = q_static_full(model.mesh.tip_dofs(3));
state.coupling_converged = converged;
state.coupling_iterations = iteration;
state.coupling_angle_error_rad = angle_error;
state.coupling_load_error = load_error;
end
