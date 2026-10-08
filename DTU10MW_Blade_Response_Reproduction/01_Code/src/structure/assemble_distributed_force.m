function force = assemble_distributed_force(model, r_load, q_flap, q_edge)

arguments
    model (1,1) struct
    r_load double
    q_flap double
    q_edge double
end

r_load = r_load(:); q_flap = q_flap(:); q_edge = q_edge(:);
if numel(r_load) ~= numel(q_flap) || numel(r_load) ~= numel(q_edge)
    error('DTU10MW:LoadSize', 'r_load, q_flap and q_edge must have equal lengths.');
end
if any(diff(r_load) <= 0)
    error('DTU10MW:LoadGrid', 'Distributed-load radius must be strictly increasing.');
end
qf = interp1(r_load,q_flap,model.mesh.mid_r,'pchip','extrap');
qe = interp1(r_load,q_edge,model.mesh.mid_r,'pchip','extrap');
F = zeros(model.mesh.ndof,1);
for e = 1:model.mesh.number_of_elements
    L = model.mesh.element_length(e);
    f_flap = qf(e)*[L/2;L^2/12;L/2;-L^2/12];
    f_edge = qe(e)*[L/2;L^2/12;L/2;-L^2/12];
    local = zeros(8,1);
    local([1 2 5 6]) = f_flap;
    local([3 4 7 8]) = f_edge;
    dof = model.mesh.element_dofs(e,:);
    F(dof) = F(dof)+local;
end
force.full = F;
force.free = F(model.mesh.free_dofs);
force.q_flap_mid = qf;
force.q_edge_mid = qe;
end
