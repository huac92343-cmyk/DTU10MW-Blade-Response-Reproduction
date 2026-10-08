function checks = validate_DTU10MW_data(data, cfg)

arguments
    data (1,1) struct
    cfg (1,1) struct
end

checks = struct();
check_radius(data.geometry.r_m, data, 'Blade_Geometry');
check_radius(data.structure.r_m, data, 'Blade_Structure');
check_radius(data.aerodynamic_stations.r_m, data, 'Aerodynamic_Stations');

positive_columns(data.geometry, {'chord_m'}, 'Blade_Geometry');
positive_columns(data.structure, {'mass_per_length_kg_m','EI_flap_Nm2', ...
    'EI_edge_Nm2','GJ_Nm2','GA_flap_N','GA_edge_N','EA_N'}, ...
    'Blade_Structure');
positive_columns(data.aerodynamic_stations, {'chord_m', ...
    'relative_thickness_pct'}, 'Aerodynamic_Stations');

if any(data.structure.elastic_axis_aft_le_chord < 0 | ...
        data.structure.elastic_axis_aft_le_chord > 1) || ...
        any(data.structure.aero_reference_axis_aft_le_chord < 0 | ...
        data.structure.aero_reference_axis_aft_le_chord > 1)
    error('DTU10MW:AxisLocation', ...
        'Axis locations aft of the leading edge must lie in [0,1] chord.');
end

mass = trapz(double(data.structure.r_m), ...
    double(data.structure.mass_per_length_kg_m));
checks.integrated_blade_mass_kg = mass;
checks.reference_blade_mass_kg = data.reference_blade_mass;
checks.mass_relative_error = abs(mass-data.reference_blade_mass) / data.reference_blade_mass;
if checks.mass_relative_error > cfg.structure.mass_relative_tolerance
    error('DTU10MW:MassMismatch', ...
        'Integrated mass %.6g kg differs from reference %.6g kg by %.3f%% (limit %.3f%%).', ...
        mass, data.reference_blade_mass, 100*checks.mass_relative_error, ...
        100*cfg.structure.mass_relative_tolerance);
end

v = double(data.operating.wind_speed_m_s);
if any(~isfinite(v)) || any(diff(v) <= 0)
    error('DTU10MW:OperatingGrid', ...
        'Operating_Curve wind_speed_m_s must be finite and strictly increasing.');
end
if data.V_rated < v(1) || data.V_rated > v(end)
    error('DTU10MW:RatedWind', 'Rated wind speed is outside the operating curve.');
end

if numel(data.station_polars) ~= height(data.aerodynamic_stations)
    error('DTU10MW:AirfoilMapping', ...
        'Each aerodynamic station must have one thickness-interpolated polar.');
end
thickness = double(data.aerodynamic_stations.relative_thickness_pct);
catalog_thickness = [data.base_polars.thickness_pct];
if any(thickness < min(catalog_thickness)-1e-10 | ...
        thickness > max(catalog_thickness)+1e-10)
    error('DTU10MW:AirfoilMapping', ...
        'Aerodynamic thickness lies outside the official polar catalog.');
end
end

function check_radius(r, data, filename)
r = double(r(:));
if any(~isfinite(r)) || any(diff(r) <= 0)
    error('DTU10MW:SpanGrid', '%s r_m must be finite and strictly increasing.', filename);
end
tol = 100*eps(max(1,data.R));
if r(1) < data.R_hub-tol || r(end) > data.R+tol
    error('DTU10MW:TipOverrun', ...
        '%s spans [%g,%g] m outside [%g,%g] m.', ...
        filename, r(1), r(end), data.R_hub, data.R);
end
end

function positive_columns(tbl, names, filename)
for k = 1:numel(names)
    x = double(tbl.(names{k}));
    if any(~isfinite(x)) || any(x <= 0)
        error('DTU10MW:NonpositiveData', '%s column %s must be finite and positive.', ...
            filename, names{k});
    end
end
end
