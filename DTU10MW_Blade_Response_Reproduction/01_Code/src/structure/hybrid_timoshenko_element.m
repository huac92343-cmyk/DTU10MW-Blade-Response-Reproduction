function element = hybrid_timoshenko_element(L, mass_per_length, ...
    EI_flap, EI_edge, GA_flap, GA_edge, axial_force, axis_angle_rad, type)

arguments
    L (1,1) double {mustBePositive}
    mass_per_length (1,1) double {mustBePositive}
    EI_flap (1,1) double {mustBePositive}
    EI_edge (1,1) double {mustBePositive}
    GA_flap (1,1) double {mustBePositive}
    GA_edge (1,1) double {mustBePositive}
    axial_force (1,1) double {mustBeNonnegative}
    axis_angle_rad (1,1) double {mustBeFinite}
    type (1,:) char
end

switch lower(type)
    case {'hybrid','hybrid_timoshenko','timoshenko'}
        phi_flap = 12*EI_flap/(GA_flap*L^2);
        phi_edge = 12*EI_edge/(GA_edge*L^2);
    case {'euler_bernoulli','eb'}
        phi_flap = 0;
        phi_edge = 0;
    otherwise
        error('DTU10MW:BeamType', 'Unknown beam element type "%s".', type);
end

k_flap = bending_stiffness(EI_flap,L,phi_flap);
k_edge = bending_stiffness(EI_edge,L,phi_edge);
k_geo = axial_geometric_stiffness(axial_force,L);
m_plane = consistent_mass(mass_per_length,L);

principal_K = zeros(8);
principal_Kg = zeros(8);
principal_M = zeros(8);
flap = [1 2 5 6];
edge = [3 4 7 8];
principal_K(flap,flap) = k_flap;
principal_K(edge,edge) = k_edge;
principal_Kg(flap,flap) = k_geo;
principal_Kg(edge,edge) = k_geo;
principal_M(flap,flap) = m_plane;
principal_M(edge,edge) = m_plane;

c = cos(axis_angle_rad);
s = sin(axis_angle_rad);
Tnode = [c 0 s 0; 0 c 0 s; -s 0 c 0; 0 -s 0 c];
T = blkdiag(Tnode,Tnode);
element.K_material = T.'*principal_K*T;
element.K_geometric = T.'*principal_Kg*T;
element.K = element.K_material+element.K_geometric;
element.M = T.'*principal_M*T;
element.shear_parameters = [phi_flap phi_edge];
element.transformation = T;
end

function k = bending_stiffness(EI,L,phi)
k = EI/(L^3*(1+phi))*[...
    12, 6*L, -12, 6*L;
    6*L, (4+phi)*L^2, -6*L, (2-phi)*L^2;
    -12, -6*L, 12, -6*L;
    6*L, (2-phi)*L^2, -6*L, (4+phi)*L^2];
end

function kg = axial_geometric_stiffness(N,L)
kg = N/(30*L)*[...
    36, 3*L, -36, 3*L;
    3*L, 4*L^2, -3*L, -L^2;
    -36, -3*L, 36, -3*L;
    3*L, -L^2, -3*L, 4*L^2];
end

function m = consistent_mass(mu,L)
m = mu*L/420*[...
    156, 22*L, 54, -13*L;
    22*L, 4*L^2, 13*L, -3*L^2;
    54, 13*L, 156, -22*L;
    -13*L, -3*L^2, -22*L, 4*L^2];
end
