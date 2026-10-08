function loads = aerodynamic_loads(rho, chord, W, phi, CL, CD, CM, cone_rad)

arguments
    rho double {mustBePositive}
    chord double {mustBePositive}
    W double {mustBeNonnegative}
    phi double
    CL double
    CD double {mustBeNonnegative}
    CM double
    cone_rad double
end

q = 0.5*rho.*W.^2;
Cn = CL.*cos(phi) + CD.*sin(phi);
Ct = CL.*sin(phi) - CD.*cos(phi);
loads.Cn = Cn;
loads.Ct = Ct;
loads.normal_load = q.*chord.*Cn.*cos(cone_rad);
loads.tangential_load = q.*chord.*Ct.*cos(cone_rad);
loads.pitching_moment_ref = q.*chord.^2.*CM;
end
