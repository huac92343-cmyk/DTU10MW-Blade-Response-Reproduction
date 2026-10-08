function damping = build_damping_matrix(cfg, data, model, modes, aero)

arguments
    cfg (1,1) struct
    data (1,1) struct
    model (1,1) struct
    modes (1,1) struct
    aero (1,1) struct
end

indices = cfg.damping.rayleigh_mode_indices(:);
if numel(indices) ~= 2 || any(indices < 1) || any(indices > numel(modes.omega_rad_s))
    error('DTU10MW:RayleighModes', 'Exactly two valid Rayleigh mode indices are required.');
end
omega = modes.omega_rad_s(indices);
zeta = data.structural_damping_ratio;
A = [1./(2*omega), omega/2];
coefficients = A\[zeta;zeta];
alpha = coefficients(1);
beta = coefficients(2);
C_structural = alpha*model.M+beta*model.K;
C = C_structural+aero.C;

modal_zeta_structural = zeros(numel(modes.omega_rad_s),1);
modal_zeta_total = zeros(numel(modes.omega_rad_s),1);
for k = 1:numel(modes.omega_rad_s)
    shape = modes.shape_free(:,k);
    denom = 2*modes.omega_rad_s(k)*(shape.'*model.M*shape);
    modal_zeta_structural(k) = (shape.'*C_structural*shape)/denom;
    modal_zeta_total(k) = (shape.'*C*shape)/denom;
end
damping.C = sparse(0.5*(C+C.'));
damping.C_structural = sparse(C_structural);
damping.C_aerodynamic = aero.C;
damping.rayleigh_alpha = alpha;
damping.rayleigh_beta = beta;
damping.target_structural_zeta = zeta;
damping.modal_zeta_structural = modal_zeta_structural;
damping.modal_zeta_total = modal_zeta_total;
end
