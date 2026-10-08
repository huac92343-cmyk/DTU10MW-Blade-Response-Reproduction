function bem = bem_solver_DTU10MW(cfg, data, wind_speed, ...
    Omega, pitch_rad, azimuth_rad, psi_elastic_rad)

arguments
    cfg (1,1) struct
    data (1,1) struct
    wind_speed (1,1) double {mustBeNonnegative}
    Omega (1,1) double {mustBeNonnegative}
    pitch_rad (1,1) double {mustBeFinite}
    azimuth_rad (1,1) double {mustBeFinite}
    psi_elastic_rad double = []
end

r = double(data.aerodynamic_stations.r_m(:));
chord = double(data.aerodynamic_stations.chord_m(:));
twist = deg2rad(double(data.aerodynamic_stations.aerodynamic_twist_deg(:)));
ids = string({data.station_polars.id}).';
n = numel(r);
radius_tolerance = 100*eps(max(data.R,1));
boundary_station = abs(r-data.R_hub) <= radius_tolerance | ...
    abs(r-data.R) <= radius_tolerance;
if isempty(psi_elastic_rad)
    psi_elastic_rad = zeros(n,1);
else
    psi_elastic_rad = psi_elastic_rad(:);
    if numel(psi_elastic_rad) ~= n
        error('DTU10MW:TwistSize', 'psi_elastic_rad must have one value per aerodynamic station.');
    end
end

kin = coordinate_transform_blade(data, azimuth_rad, [-wind_speed;0;0]);
Vfree_flap = -kin.wind_blade(1);
Vfree_edge = -kin.wind_blade(2);
if Vfree_flap < 0
    error('DTU10MW:WindDirection', ...
        'Configured wind/coordinate convention gives negative axial inflow %.6g m/s.', Vfree_flap);
end

[a, a_prime] = deal(zeros(n,1));
[iterations, converged] = deal(zeros(n,1), false(n,1));
for i = 1:n
    if boundary_station(i)
        a(i) = 0;
        a_prime(i) = 0;
        converged(i) = true;
        iterations(i) = 0;
        continue
    end
    sigma = data.B*chord(i)/(2*pi*r(i));
    ai = 0.30;
    api = 0.0;
    for it = 1:cfg.aero.max_iterations
        state = section_state(ai, api, i);
        F = prandtl_loss(data.B, r(i), data.R, data.R_hub, ...
            state.phi, cfg.aero.minimum_loss_factor);
        s = sin(state.phi); c = cos(state.phi);
        denom_ap = 4*F*s*c - sigma*state.Ct;
        if abs(denom_ap) < eps
            error('DTU10MW:BEMDenominator', 'Singular tangential induction equation at r=%g m.', r(i));
        end
        ap_new = sigma*state.Ct/denom_ap;

        CT = sigma*(1-ai)^2*state.Cn/max(s^2, eps);
        if CT > 0.96*F
            disc = CT*(50-36*F) + 12*F*(3*F-4);
            if disc < -cfg.aero.discriminant_tolerance
                error('DTU10MW:HighInductionDiscriminant', ...
                    'Negative high-induction discriminant %.6g at r=%g m.', disc, r(i));
            end
            disc = max(disc,0);
            denom_hi = 36*F-50;
            if abs(denom_hi) < eps
                error('DTU10MW:HighInductionDenominator', ...
                    'Singular high-induction denominator at r=%g m.', r(i));
            end
            a_new = (18*F-20-3*sqrt(disc))/denom_hi;
        else
            momentum_disc = 1-CT/F;
            if momentum_disc < -cfg.aero.discriminant_tolerance
                error('DTU10MW:MomentumDiscriminant', ...
                    'Negative momentum discriminant %.6g at r=%g m.', ...
                    momentum_disc, r(i));
            end
            a_new = 0.5*(1-sqrt(max(momentum_disc,0)));
        end

        if ~isfinite(a_new) || ~isfinite(ap_new)
            error('DTU10MW:BEMNonfinite', 'Nonfinite induction factor at r=%g m.', r(i));
        end
        a_new = min(max(a_new,-0.2),0.95);
        ap_new = min(max(ap_new,-0.95),0.95);
        au = ai + cfg.aero.relaxation*(a_new-ai);
        apu = api + cfg.aero.relaxation*(ap_new-api);
        if max(abs([au-ai, apu-api])) < cfg.aero.tolerance
            ai = au;
            api = apu;
            converged(i) = true;
            iterations(i) = it;
            break
        end
        ai = au;
        api = apu;
        iterations(i) = it;
    end
    a(i) = ai;
    a_prime(i) = api;
end

if cfg.aero.fail_on_nonconvergence && ~all(converged)
    bad = find(~converged,1);
    error('DTU10MW:BEMNonconvergence', ...
        'BEM failed at r=%g m after %d iterations.', r(bad), iterations(bad));
end

fields = {'Vax','Vrot','phi','alpha','CL','CD','CM','dCL_dalpha', ...
    'dCD_dalpha','dCM_dalpha','Cn','Ct','W','normal_load', ...
    'tangential_load','pitching_moment_ref'};
for k = 1:numel(fields)
    bem.(fields{k}) = zeros(n,1);
end
for i = 1:n
    final = section_state(a(i), a_prime(i), i);
    if boundary_station(i)
        final.normal_load = 0;
        final.tangential_load = 0;
        final.pitching_moment_ref = 0;
    end
    names = fieldnames(final);
    for k = 1:numel(names)
        if isfield(bem,names{k})
            bem.(names{k})(i) = final.(names{k});
        end
    end
end

bem.r = r;
bem.chord = chord;
bem.structural_twist_rad = twist;
bem.psi_elastic_rad = psi_elastic_rad;
bem.effective_angle_rad = twist + pitch_rad + psi_elastic_rad;
bem.airfoil_id = ids;
bem.a = a;
bem.a_prime = a_prime;
bem.bem_converged = converged;
bem.bem_iterations = iterations;
bem.coordinate_state = kin;
bem.boundary_zero_load = boundary_station;

    function state = section_state(ai, api, idx)
        state.Vax = Vfree_flap*(1-ai);
        state.Vrot = Omega*r(idx)*cos(data.cone_rad)*(1+api) - Vfree_edge;
        state.phi = atan2(state.Vax, state.Vrot);
        state.alpha = state.phi - (twist(idx)+pitch_rad+psi_elastic_rad(idx));
        coeff = interpolate_airfoil_coefficients(data.station_polars(idx), ...
            state.alpha, cfg);
        state.CL = coeff.CL;
        state.CD = coeff.CD;
        state.CM = coeff.CM;
        state.dCL_dalpha = coeff.dCL_dalpha;
        state.dCD_dalpha = coeff.dCD_dalpha;
        state.dCM_dalpha = coeff.dCM_dalpha;
        state.W = hypot(state.Vax,state.Vrot);
        ld = aerodynamic_loads(data.rho, chord(idx), state.W, state.phi, ...
            state.CL, state.CD, state.CM, data.cone_rad);
        state.Cn = ld.Cn;
        state.Ct = ld.Ct;
        state.normal_load = ld.normal_load;
        state.tangential_load = ld.tangential_load;
        state.pitching_moment_ref = ld.pitching_moment_ref;
    end
end

function F = prandtl_loss(B, r, R, R_hub, phi, minimum_F)
if r < R_hub || r > R
    error('DTU10MW:TipOverrun', 'BEM radius %.9g m is outside [%g,%g] m.', r, R_hub, R);
end
s = max(abs(sin(phi)),1e-8);
r_eff = min(max(r,R_hub+eps(R_hub)),R-eps(R));
f_tip = B*(R-r_eff)/(2*r_eff*s);
f_hub = B*(r_eff-R_hub)/(2*R_hub*s);
F_tip = (2/pi)*acos(min(max(exp(-f_tip),0),1));
F_hub = (2/pi)*acos(min(max(exp(-f_hub),0),1));
F = max(F_tip*F_hub, minimum_F);
end
