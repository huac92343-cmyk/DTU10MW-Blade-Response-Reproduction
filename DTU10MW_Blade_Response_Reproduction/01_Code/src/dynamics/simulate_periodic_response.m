function response = simulate_periodic_response(cfg, periodic)

arguments
    cfg (1,1) struct
    periodic (1,1) struct
end

period = 2*pi/periodic.Omega_rad_s;
steps_per_revolution = max(12,ceil(period/cfg.time.dt));
dt = period/steps_per_revolution;
maximum_steps = cfg.time.maximum_revolutions*steps_per_revolution;
time = (0:maximum_steps)*dt;
force_fun = @(t) periodic_force_at_time(t,periodic.Omega_rad_s, ...
    periodic.force_samples);
q0 = periodic.model.K\periodic.mean_force;
v0 = zeros(size(q0));
history = houbolt_integrate(periodic.model.M,periodic.damping.C, ...
    periodic.model.K,force_fun,time,q0,v0);
if any(~isfinite(history.displacement),'all') || ...
        any(~isfinite(history.velocity),'all') || ...
        any(~isfinite(history.acceleration),'all')
    error('DTU10MW:NonfiniteResponse', ...
        'The time integration produced nonfinite displacement, velocity, or acceleration.');
end

tip_flap_full = periodic.model.mesh.tip_dofs(1);
tip_edge_full = periodic.model.mesh.tip_dofs(3);
tip_flap_free = find(periodic.model.mesh.free_dofs == tip_flap_full,1);
tip_edge_free = find(periodic.model.mesh.free_dofs == tip_edge_full,1);
tip_all = history.displacement([tip_flap_free tip_edge_free],:);
converged_revolution = cfg.time.maximum_revolutions;
periodic_error = inf;
for revolution = cfg.time.minimum_revolutions:cfg.time.maximum_revolutions
    current = (revolution-1)*steps_per_revolution+1: ...
        revolution*steps_per_revolution+1;
    previous = current-steps_per_revolution;
    periodic_error = norm(tip_all(:,current)-tip_all(:,previous),'fro')/ ...
        max(norm(tip_all(:,previous),'fro'),1e-8);
    if periodic_error <= cfg.time.periodic_relative_tolerance
        converged_revolution = revolution;
        break
    end
end
if ~isfinite(periodic_error) || periodic_error > cfg.time.periodic_relative_tolerance
    error('DTU10MW:PeriodicNonconvergence', ...
        'Periodic response failed after %d revolutions: relative error %.6g exceeds %.6g.', ...
        cfg.time.maximum_revolutions,periodic_error,cfg.time.periodic_relative_tolerance);
end
last_index = converged_revolution*steps_per_revolution+1;
history.time_s = history.time_s(1:last_index);
history.displacement = history.displacement(:,1:last_index);
history.velocity = history.velocity(:,1:last_index);
history.acceleration = history.acceleration(:,1:last_index);

steady_revolutions = min(cfg.time.steady_cycles,converged_revolution);
first_steady = (converged_revolution-steady_revolutions)* ...
    steps_per_revolution+1;
steady_index = first_steady:last_index;
tip_flap = history.displacement(tip_flap_free,:);
tip_edge = history.displacement(tip_edge_free,:);

response.history = history;
response.tip_flap_m = tip_flap;
response.tip_edge_m = tip_edge;
response.steady_index = steady_index;
response.steps_per_revolution = steps_per_revolution;
response.revolutions_completed = converged_revolution;
response.periodic_relative_error = periodic_error;
response.periodic_converged = periodic_error <= cfg.time.periodic_relative_tolerance;
response.metrics_flap = response_metrics(tip_flap(steady_index));
response.metrics_edge = response_metrics(tip_edge(steady_index));
end

function metrics = response_metrics(x)
metrics.mean = mean(x);
metrics.maximum = max(x);
metrics.minimum = min(x);
metrics.peak_absolute = max(abs(x));
metrics.amplitude = 0.5*(max(x)-min(x));
metrics.rms = rms(x);
end
