function coeff = interpolate_airfoil_coefficients(polar, alpha_rad, cfg)

arguments
    polar (1,1) struct
    alpha_rad double
    cfg (1,1) struct
end

alpha_deg = rad2deg(alpha_rad);
amin = polar.alpha_deg(1);
amax = polar.alpha_deg(end);
outside = alpha_deg < amin | alpha_deg > amax;
if any(outside, 'all')
    switch lower(cfg.aero.polar_out_of_range_policy)
        case 'error'
            bad = alpha_deg(find(outside,1));
            error('DTU10MW:PolarRange', ...
                'Angle %.6g deg is outside [%g,%g] deg for airfoil %s; unconstrained extrapolation is disabled.', ...
                bad, amin, amax, polar.id);
        otherwise
            error('DTU10MW:PolarPolicy', ...
                'Unsupported polar_out_of_range_policy "%s".', ...
                cfg.aero.polar_out_of_range_policy);
    end
end

coeff.CL = interp1(polar.alpha_deg, polar.CL, alpha_deg, 'pchip');
coeff.CD = interp1(polar.alpha_deg, polar.CD, alpha_deg, 'pchip');
coeff.CM = interp1(polar.alpha_deg, polar.CM, alpha_deg, 'pchip');

coeff.dCL_dalpha = local_slope(polar.alpha_deg, polar.CL, alpha_deg) * (180/pi);
coeff.dCD_dalpha = local_slope(polar.alpha_deg, polar.CD, alpha_deg) * (180/pi);
coeff.dCM_dalpha = local_slope(polar.alpha_deg, polar.CM, alpha_deg) * (180/pi);
end

function dydx = local_slope(x, y, xq)
d = gradient(y(:), x(:));
dydx = interp1(x, d, xq, 'linear');
end
