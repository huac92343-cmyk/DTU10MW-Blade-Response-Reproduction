function model = assemble_blade_matrices(cfg, data, Omega, pitch_rad, ...
    elastic_twist_rad)

arguments
    cfg (1,1) struct
    data (1,1) struct
    Omega (1,1) double {mustBeNonnegative}
    pitch_rad (1,1) double {mustBeFinite}
    elastic_twist_rad double = []
end

mesh = build_blade_mesh(data,cfg.structure.number_of_elements);
r_node = mesh.r;
r_mid = mesh.mid_r;
s = data.structure;
mass_node = interp1(s.r_m,s.mass_per_length_kg_m,r_node,'pchip');
centrifugal_density = mass_node.*Omega.^2.*r_node*cos(data.cone_rad);
axial_node = zeros(size(r_node));
for j = numel(r_node)-1:-1:1
    axial_node(j) = axial_node(j+1)+0.5*(centrifugal_density(j)+ ...
        centrifugal_density(j+1))*(r_node(j+1)-r_node(j));
end

mass = interp1(s.r_m,s.mass_per_length_kg_m,r_mid,'pchip');
GA_flap = interp1(s.r_m,s.GA_flap_N,r_mid,'pchip');
GA_edge = interp1(s.r_m,s.GA_edge_N,r_mid,'pchip');

if isempty(elastic_twist_rad)
    elastic_mid = zeros(size(r_mid));
elseif isscalar(elastic_twist_rad)
    elastic_mid = repmat(elastic_twist_rad,size(r_mid));
elseif numel(elastic_twist_rad) == height(s)
    elastic_mid = interp1(s.r_m,elastic_twist_rad(:),r_mid,'pchip');
elseif numel(elastic_twist_rad) == mesh.number_of_nodes
    elastic_mid = 0.5*(elastic_twist_rad(1:end-1)+elastic_twist_rad(2:end));
else
    error('DTU10MW:TwistSize', ...
        'Elastic twist must be scalar, structural-grid, or mesh-node sized.');
end
section_stiffness = transform_section_stiffness(data,r_mid,pitch_rad,elastic_mid);
EI_flap = section_stiffness.EI_principal_x_Nm2;
EI_edge = section_stiffness.EI_principal_y_Nm2;
axis_angle = section_stiffness.total_axis_angle_rad;
axial_mid = 0.5*(axial_node(1:end-1)+axial_node(2:end));

K_material = zeros(mesh.ndof);
K_geometric = zeros(mesh.ndof);
M = zeros(mesh.ndof);
shear_parameters = zeros(mesh.number_of_elements,2);
for e = 1:mesh.number_of_elements
    elem = hybrid_timoshenko_element(mesh.element_length(e),mass(e), ...
        EI_flap(e),EI_edge(e),GA_flap(e),GA_edge(e),axial_mid(e), ...
        axis_angle(e),cfg.structure.element_type);
    dof = mesh.element_dofs(e,:);
    K_material(dof,dof) = K_material(dof,dof)+elem.K_material;
    K_geometric(dof,dof) = K_geometric(dof,dof)+elem.K_geometric;
    M(dof,dof) = M(dof,dof)+elem.M;
    shear_parameters(e,:) = elem.shear_parameters;
end
K = K_material+K_geometric;
free = mesh.free_dofs;
symmetry_error_K = norm(K-K.','fro')/max(norm(K,'fro'),eps);
symmetry_error_M = norm(M-M.','fro')/max(norm(M,'fro'),eps);
if max(symmetry_error_K,symmetry_error_M) > cfg.structure.symmetry_tolerance
    error('DTU10MW:MatrixSymmetry', ...
        'Beam matrix symmetry error K=%.3g, M=%.3g.', ...
        symmetry_error_K,symmetry_error_M);
end
[~,mass_chol_flag] = chol(M(free,free));
[~,stiffness_chol_flag] = chol(K(free,free));
if mass_chol_flag ~= 0 || stiffness_chol_flag ~= 0
    error('DTU10MW:MatrixPositiveDefinite', ...
        'Constrained beam matrices are not positive definite (M=%d,K=%d).', ...
        mass_chol_flag,stiffness_chol_flag);
end

model.mesh = mesh;
model.K_full = sparse(K);
model.K_material_full = sparse(K_material);
model.K_geometric_full = sparse(K_geometric);
model.M_full = sparse(M);
model.K = model.K_full(free,free);
model.K_material = model.K_material_full(free,free);
model.K_geometric = model.K_geometric_full(free,free);
model.M = model.M_full(free,free);
model.axial_force_node_N = axial_node;
model.axis_angle_mid_rad = axis_angle;
model.shear_parameters = shear_parameters;
model.Omega = Omega;
model.pitch_rad = pitch_rad;
model.elastic_twist_mid_rad = elastic_mid;
model.section_stiffness = section_stiffness;
model.symmetry_error_K = symmetry_error_K;
model.symmetry_error_M = symmetry_error_M;
end
