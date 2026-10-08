function radius = compute_equivalent_radius(data, periodic, response)

arguments
    data (1,1) struct
    periodic (1,1) struct
    response (1,1) struct
end

b = data.cone_rad;
prebend_flap = interp1(data.geometry.r_m,data.geometry.prebend_flap_m, ...
    data.R,'linear');
prebend_edge = interp1(data.geometry.r_m,data.geometry.prebend_edge_m, ...
    data.R,'linear');
u_flap = response.tip_flap_m;
u_edge = response.tip_edge_m;
radial_component = data.R*cos(b)-(prebend_flap+u_flap)*sin(b);
edge_component = prebend_edge+u_edge;
dynamic = hypot(radial_component,edge_component);
static = hypot(data.R*cos(b)-prebend_flap*sin(b),prebend_edge);

radius.nominal_m = data.R;
radius.coned_undeformed_m = data.R*cos(b);
radius.static_geometry_m = static;
radius.dynamic_m = dynamic;
radius.steady_dynamic_m = dynamic(response.steady_index);
radius.mean_m = mean(radius.steady_dynamic_m);
radius.minimum_m = min(radius.steady_dynamic_m);
radius.maximum_m = max(radius.steady_dynamic_m);
radius.periodic_azimuth_samples = periodic.azimuth_rad;
end
