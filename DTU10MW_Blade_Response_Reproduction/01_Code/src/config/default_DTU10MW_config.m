function cfg = default_DTU10MW_config()
% Numerical settings are defined here. The Excel Code_Parameters sheet
% lists these settings for inspection; it does not override this configuration.

source_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
package_root = fileparts(source_root);

cfg.paths.root_dir = package_root;
cfg.paths.data_dir = fullfile(package_root, '02_Input_Parameters');
cfg.paths.results_dir = fullfile(package_root, '03_Results', 'MATLAB_Data');
cfg.paths.figure_dir = fullfile(package_root, '03_Results', 'Figures');

cfg.coordinates.blade_axis = '+r: shaft centre toward blade tip';
cfg.coordinates.flap_direction = '+x_b: nominal downwind/flapwise direction';
cfg.coordinates.edge_direction = '+y_b: direction of rotation at zero deformation';
cfg.coordinates.wake_x = '+x_w: downstream along the nominal shaft axis';
cfg.coordinates.wake_y = '+y_w: horizontal transverse direction';
cfg.coordinates.wake_z = '+z_w: vertically upward';
cfg.coordinates.positive_twist = 'nose toward positive pitch about +r';
cfg.coordinates.positive_pitch = 'feathering rotation about +r';
cfg.coordinates.azimuth_zero = 'blade vertically upward';
cfg.coordinates.azimuth_positive = 'direction of rotor rotation';

cfg.structure.element_type = 'hybrid';
cfg.structure.number_of_elements = 60;
cfg.structure.mode_count = 6;
cfg.structure.mass_relative_tolerance = 0.02;
cfg.structure.symmetry_tolerance = 1e-10;

cfg.aero.max_iterations = 300;
cfg.aero.tolerance = 1e-7;
cfg.aero.relaxation = 0.30;
cfg.aero.discriminant_tolerance = 1e-10;
cfg.aero.minimum_loss_factor = 1e-4;
cfg.aero.fail_on_nonconvergence = true;
cfg.aero.polar_out_of_range_policy = 'error';

cfg.coupling.max_iterations = 50;
cfg.coupling.tolerance_rad = 1e-5;
cfg.coupling.load_relative_tolerance = 1e-4;
cfg.coupling.relaxation = 0.35;
cfg.coupling.fail_on_nonconvergence = true;
cfg.twist_coupling_factor = 4.0;

cfg.damping.rayleigh_mode_indices = [1 3];

cfg.time.dt = 0.04;
cfg.time.azimuth_samples = 12;
cfg.time.minimum_revolutions = 20;
cfg.time.maximum_revolutions = 40;
cfg.time.steady_cycles = 3;
cfg.time.periodic_relative_tolerance = 0.01;

cfg.output.save_figures = true;
cfg.output.figure_formats = {'png','pdf','fig'};

cfg.analysis.wind_speeds_m_s = [];

cfg.validation.mesh_elements = [30 60 90 120];
cfg.validation.time_step_factors = [1 0.5 0.25];

cfg.software.matlab_release = 'R2023b';
cfg.software.code_version = 'DTU10MW-model-1.0.0';
end
