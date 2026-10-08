function data = load_DTU10MW_data(cfg)

arguments
    cfg (1,1) struct
end

data_dir = cfg.paths.data_dir;
required_files = {
    '01_Turbine_Operating_Parameters.xlsx'
    '02_Blade_Geometry_Aerodynamics.xlsx'
    '03_Blade_Structural_Parameters.xlsx'
    '04_Airfoil_Polar_Data.xlsx'};

missing = required_files(~cellfun(@(f) isfile(fullfile(data_dir,f)), required_files));
polar_dir = fullfile(data_dir, '04_Airfoil_Polar_Data.xlsx');
if ~isempty(missing)
    error('DTU10MW:MissingData', ...
        'Missing DTU 10 MW input(s) under %s: %s. Restore the packaged runtime data files.', ...
        data_dir, strjoin(missing, ', '));
end

manifest_cells = readcell(fullfile(data_dir, ...
    '01_Turbine_Operating_Parameters.xlsx'), 'Sheet','Turbine');
manifest = struct();
for k = 2:size(manifest_cells,1)
    manifest.(char(string(manifest_cells{k,1}))) = manifest_cells{k,2};
end
required_manifest = {'model_name','R_m','R_hub_m','blade_count', ...
    'rated_wind_speed_m_s','cone_deg','shaft_tilt_deg','yaw_deg', ...
    'air_density_kg_m3','blade_reference_mass_kg', ...
    'blade_structural_log_decrement_percent'};
assert_fields(manifest, required_manifest, 'Turbine');
if ~contains(upper(string(manifest.model_name)), 'DTU') || ...
        ~contains(upper(string(manifest.model_name)), '10')
    error('DTU10MW:WrongModel', ...
        'Turbine model_name must identify a DTU 10 MW model, got "%s".', ...
        string(manifest.model_name));
end

geom = readtable(fullfile(data_dir, '02_Blade_Geometry_Aerodynamics.xlsx'), ...
    'Sheet','Blade_Geometry','VariableNamingRule','preserve');
structure = readtable(fullfile(data_dir, '03_Blade_Structural_Parameters.xlsx'), ...
    'Sheet','Blade_Structure','VariableNamingRule','preserve');
operating = readtable(fullfile(data_dir, '01_Turbine_Operating_Parameters.xlsx'), ...
    'Sheet','Operating_Curve','VariableNamingRule','preserve');
aero_stations = readtable(fullfile(data_dir, '02_Blade_Geometry_Aerodynamics.xlsx'), ...
    'Sheet','Aerodynamic_Stations','VariableNamingRule','preserve');
polar_catalog = readtable(polar_dir, 'Sheet','Polar_Catalog', ...
    'VariableNamingRule','preserve','TextType','string');

require_columns(geom, {'r_m','chord_m','aerodynamic_twist_deg', ...
    'prebend_flap_m','prebend_edge_m','relative_thickness_pct', ...
    'pitch_axis_aft_le_chord'}, 'Blade_Geometry');
require_columns(structure, {'r_m','mass_per_length_kg_m','EI_flap_Nm2', ...
    'EI_edge_Nm2','GJ_Nm2','GA_flap_N','GA_edge_N','EA_N', ...
    'structural_pitch_deg','elastic_axis_aft_le_chord', ...
    'aero_reference_axis_aft_le_chord'}, ...
    'Blade_Structure');
require_columns(operating, {'wind_speed_m_s','rpm','pitch_deg'}, ...
    'Operating_Curve');
require_columns(aero_stations, {'r_m','chord_m','aerodynamic_twist_deg', ...
    'relative_thickness_pct','polar_set_id','pitch_axis_aft_le_chord'}, ...
    'Aerodynamic_Stations');
require_columns(polar_catalog, {'airfoil_id','relative_thickness_pct', ...
    'polar_sheet'}, 'Polar_Catalog');

data.manifest = manifest;
data.geometry = geom;
data.structure = structure;
data.operating = operating;
data.aerodynamic_stations = aero_stations;
data.polar_catalog = polar_catalog;
[data.base_polars, data.station_polars] = load_airfoil_polars( ...
    polar_dir, polar_catalog, aero_stations);

data.R = double(manifest.R_m);
data.R_hub = double(manifest.R_hub_m);
data.B = double(manifest.blade_count);
data.rho = double(manifest.air_density_kg_m3);
data.g = 9.80665;
data.V_rated = double(manifest.rated_wind_speed_m_s);
data.cone_rad = deg2rad(double(manifest.cone_deg));
data.shaft_tilt_rad = deg2rad(double(manifest.shaft_tilt_deg));
data.yaw_rad = deg2rad(double(manifest.yaw_deg));
data.reference_blade_mass = double(manifest.blade_reference_mass_kg);
delta = double(manifest.blade_structural_log_decrement_percent)/100;
data.structural_damping_ratio = delta/sqrt((2*pi)^2+delta^2);

data.validation = validate_DTU10MW_data(data, cfg);
end

function require_columns(tbl, names, filename)
missing = names(~ismember(names, tbl.Properties.VariableNames));
if ~isempty(missing)
    error('DTU10MW:MissingColumns', '%s lacks column(s): %s.', ...
        filename, strjoin(missing, ', '));
end
end

function assert_fields(s, names, filename)
missing = names(~isfield(s,names));
if ~isempty(missing)
    error('DTU10MW:MissingFields', '%s lacks field(s): %s.', ...
        filename, strjoin(missing, ', '));
end
end
