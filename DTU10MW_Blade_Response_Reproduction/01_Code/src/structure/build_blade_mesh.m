function mesh = build_blade_mesh(data, number_of_elements)

arguments
    data (1,1) struct
    number_of_elements (1,1) double {mustBeInteger,mustBePositive}
end

mesh.number_of_elements = number_of_elements;
mesh.number_of_nodes = number_of_elements+1;
mesh.r = linspace(data.R_hub, data.R, mesh.number_of_nodes).';
mesh.element_length = diff(mesh.r);
mesh.mid_r = 0.5*(mesh.r(1:end-1)+mesh.r(2:end));
mesh.dofs_per_node = 4;
mesh.ndof = mesh.dofs_per_node*mesh.number_of_nodes;
mesh.fixed_dofs = (1:4).';
mesh.free_dofs = (5:mesh.ndof).';
mesh.tip_dofs = mesh.ndof-3:mesh.ndof;

mesh.element_dofs = zeros(number_of_elements,8);
for e = 1:number_of_elements
    mesh.element_dofs(e,:) = [4*e-3:4*e, 4*e+1:4*e+4];
end
end
