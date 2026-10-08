function damping = assemble_aerodynamic_damping(data, model, bem)

arguments
    data (1,1) struct
    model (1,1) struct
    bem (1,1) struct
end

phi = bem.phi;
W = max(bem.W,1e-8);
CL = bem.CL; CD = bem.CD;
dCL = bem.dCL_dalpha; dCD = bem.dCD_dalpha;
Cn = bem.Cn; Ct = bem.Ct;
cphi = cos(phi); sphi = sin(phi);
dCn_dphi = dCL.*cphi-CL.*sphi+dCD.*sphi+CD.*cphi;
dCt_dphi = dCL.*sphi+CL.*cphi-dCD.*cphi+CD.*sphi;
factor = data.rho.*bem.chord.*cos(data.cone_rad);
dN_dVax = factor.*W.*(sphi.*Cn+0.5*cphi.*dCn_dphi);
dT_dVrot = factor.*W.*(cphi.*Ct-0.5*sphi.*dCt_dphi);

c_flap = dN_dVax;
c_edge = -dT_dVrot;
if isfield(bem,'boundary_zero_load')
    c_flap(bem.boundary_zero_load) = 0;
    c_edge(bem.boundary_zero_load) = 0;
end
c_flap_mid = interp1(bem.r,c_flap,model.mesh.mid_r,'pchip');
c_edge_mid = interp1(bem.r,c_edge,model.mesh.mid_r,'pchip');
C = zeros(model.mesh.ndof);
for e = 1:model.mesh.number_of_elements
    L = model.mesh.element_length(e);
    cf = c_flap_mid(e)*L/6*[2 1;1 2];
    ce = c_edge_mid(e)*L/6*[2 1;1 2];
    dof = model.mesh.element_dofs(e,:);
    C(dof([1 5]),dof([1 5])) = C(dof([1 5]),dof([1 5]))+cf;
    C(dof([3 7]),dof([3 7])) = C(dof([3 7]),dof([3 7]))+ce;
end
damping.C_full = sparse(0.5*(C+C.'));
damping.C = damping.C_full(model.mesh.free_dofs,model.mesh.free_dofs);
damping.c_flap_Ns_m2 = c_flap;
damping.c_edge_Ns_m2 = c_edge;
damping.c_flap_mid_Ns_m2 = c_flap_mid;
damping.c_edge_mid_Ns_m2 = c_edge_mid;
damping.minimum_flap = min(c_flap);
damping.minimum_edge = min(c_edge);
end
