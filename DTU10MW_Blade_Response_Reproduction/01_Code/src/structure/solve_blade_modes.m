function modes = solve_blade_modes(model, mode_count)

arguments
    model (1,1) struct
    mode_count (1,1) double {mustBeInteger,mustBePositive}
end

[V,D] = eig(full(model.K),full(model.M),'vector');
valid = isfinite(D) & real(D) > 0 & abs(imag(D)) < 1e-8*max(1,abs(real(D)));
D = real(D(valid));
V = real(V(:,valid));
[D,order] = sort(D,'ascend');
V = V(:,order);
n = min(mode_count,numel(D));
D = D(1:n);
V = V(:,1:n);

for k = 1:n
    V(:,k) = V(:,k)/sqrt(V(:,k).'*model.M*V(:,k));
end
full_shapes = zeros(model.mesh.ndof,n);
full_shapes(model.mesh.free_dofs,:) = V;
flap_dof = 1:4:model.mesh.ndof;
edge_dof = 3:4:model.mesh.ndof;
flap_fraction = zeros(n,1);
classification = strings(n,1);
for k = 1:n
    ef = sum(full_shapes(flap_dof,k).^2);
    ee = sum(full_shapes(edge_dof,k).^2);
    flap_fraction(k) = ef/max(ef+ee,eps);
    if flap_fraction(k) >= 0.5
        classification(k) = "flap-dominant";
    else
        classification(k) = "edge-dominant";
    end
end
modes.frequency_hz = sqrt(D)/(2*pi);
modes.omega_rad_s = sqrt(D);
modes.shape_free = V;
modes.shape_full = full_shapes;
modes.flap_fraction = flap_fraction;
modes.classification = classification;
end
