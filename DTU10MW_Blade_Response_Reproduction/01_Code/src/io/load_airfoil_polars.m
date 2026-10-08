function [base_polars, station_polars] = load_airfoil_polars( ...
    polar_dir, catalog, aero_stations)

arguments
    polar_dir (1,:) char
    catalog table
    aero_stations table
end

if numel(unique(string(catalog.airfoil_id))) ~= height(catalog) || ...
        numel(unique(string(catalog.polar_sheet))) ~= height(catalog)
    error('DTU10MW:AirfoilMapping', ...
        'Polar_Catalog must contain unique airfoil_id and polar_sheet values.');
end

base_polars = struct('id',{},'thickness_pct',{},'source_file',{}, ...
    'source_description',{},'alpha_deg',{},'CL',{},'CD',{},'CM',{});
for k = 1:height(catalog)
    id = string(catalog.airfoil_id(k));
    rel = string(catalog.polar_sheet(k));
    file = polar_dir;
    if ~isfile(file)
        error('DTU10MW:MissingPolar', 'Missing polar file for %s: %s.', id, file);
    end
    t = readtable(file, 'Sheet',char(rel),'VariableNamingRule','preserve');
    required = {'alpha_deg','CL','CD','CM'};
    missing = required(~ismember(required,t.Properties.VariableNames));
    if ~isempty(missing)
        error('DTU10MW:MissingCM', ...
            'Polar %s lacks required column(s): %s. CM may not default to zero.', ...
            file, strjoin(missing, ', '));
    end
    alpha = double(t.alpha_deg(:));
    [alpha, order] = sort(alpha);
    if any(~isfinite(alpha)) || any(diff(alpha) <= 0)
        error('DTU10MW:PolarGrid', 'alpha_deg in %s must be finite and strictly increasing.', file);
    end
    p.id = char(id);
    p.thickness_pct = double(catalog.relative_thickness_pct(k));
    p.source_file = ['04_Airfoil_Polar_Data.xlsx#' char(rel)];
    p.source_description = '';
    p.alpha_deg = alpha;
    p.CL = double(t.CL(order));
    p.CD = double(t.CD(order));
    p.CM = double(t.CM(order));
    if any(~isfinite([p.CL; p.CD; p.CM]))
        error('DTU10MW:PolarValues', 'Polar coefficients in %s must be finite.', file);
    end
    base_polars(end+1) = p; %#ok<AGROW>
end

[~,order] = sort([base_polars.thickness_pct]);
base_polars = base_polars(order);
thickness_grid = [base_polars.thickness_pct].';
alpha = base_polars(1).alpha_deg;
for k = 2:numel(base_polars)
    if numel(base_polars(k).alpha_deg) ~= numel(alpha) || ...
            max(abs(base_polars(k).alpha_deg-alpha)) > 1e-10
        error('DTU10MW:PolarGrid', ...
            'Official thickness profiles must share a common alpha grid.');
    end
end

CL = [base_polars.CL];
CD = [base_polars.CD];
CM = [base_polars.CM];
t_station = double(aero_stations.relative_thickness_pct(:));
if any(t_station < thickness_grid(1)-1e-10 | ...
        t_station > thickness_grid(end)+1e-10)
    error('DTU10MW:AirfoilThickness', ...
        'Aerodynamic-station thickness is outside the official polar range.');
end

station_polars = repmat(struct('id','','thickness_pct',0, ...
    'lower_profile','','upper_profile','','blend_weight',0, ...
    'alpha_deg',alpha,'CL',alpha,'CD',alpha,'CM',alpha), ...
    numel(t_station),1);
for i = 1:numel(t_station)
    t = t_station(i);
    upper = find(thickness_grid >= t-1e-12,1,'first');
    lower = find(thickness_grid <= t+1e-12,1,'last');
    if lower == upper
        weight = 0;
    else
        weight = (t-thickness_grid(lower))/ ...
            (thickness_grid(upper)-thickness_grid(lower));
    end
    station_polars(i).id = sprintf('station_%03d',i);
    station_polars(i).thickness_pct = t;
    station_polars(i).lower_profile = base_polars(lower).id;
    station_polars(i).upper_profile = base_polars(upper).id;
    station_polars(i).blend_weight = weight;
    station_polars(i).alpha_deg = alpha;
    station_polars(i).CL = (1-weight)*CL(:,lower)+weight*CL(:,upper);
    station_polars(i).CD = (1-weight)*CD(:,lower)+weight*CD(:,upper);
    station_polars(i).CM = (1-weight)*CM(:,lower)+weight*CM(:,upper);
end
end
