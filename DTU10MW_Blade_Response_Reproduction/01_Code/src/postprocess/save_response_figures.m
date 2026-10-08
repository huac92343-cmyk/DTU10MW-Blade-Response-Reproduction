function files = save_response_figures(cfg, summary, modal_summary, ...
    rated_stiffness)

arguments
    cfg (1,1) struct
    summary table
    modal_summary table
    rated_stiffness table
end

assert(istable(modal_summary));

if ~isfolder(cfg.paths.figure_dir), mkdir(cfg.paths.figure_dir); end
files = strings(0,1);

f = figure('Visible','off','Color','w','Position',[100 100 900 620]);
plot(summary.wind_speed_m_s,summary.flap_mean_m,'o-','LineWidth',1.4);
grid on;
xlabel('Wind speed (m/s)');
ylabel('Flapwise tip displacement (m)');
title('Mean flapwise tip displacement');
files = [files; export_figure(f,cfg,'01_flapwise_tip_displacement')];
close(f);

f = figure('Visible','off','Color','w','Position',[100 100 900 620]);
plot(summary.wind_speed_m_s,summary.edge_mean_m,'o-','LineWidth',1.4);
grid on;
xlabel('Wind speed (m/s)');
ylabel('Edgewise tip displacement (m)');
title('Mean edgewise tip displacement');
files = [files; export_figure(f,cfg,'02_edgewise_tip_displacement')];
close(f);

f = figure('Visible','off','Color','w','Position',[100 100 900 620]);
plot(summary.wind_speed_m_s,summary.effective_radius_mean_m, ...
    'o-','LineWidth',1.4);
grid on;
xlabel('Wind speed (m/s)');
ylabel('Equivalent rotor radius (m)');
title('Mean equivalent rotor radius');
files = [files; export_figure(f,cfg,'03_equivalent_rotor_radius')];
close(f);

f = figure('Visible','off','Color','w','Position',[100 100 900 620]);
plot(summary.wind_speed_m_s,summary.tip_twist_total_deg, ...
    'o-','LineWidth',1.4);
grid on;
xlabel('Wind speed (m/s)');
ylabel('Tip twist angle (deg)');
title('Total elastic tip twist');
files = [files; export_figure(f,cfg,'04_tip_twist_angle')];
close(f);

if ~isempty(rated_stiffness)
    f = figure('Visible','off','Color','w','Position',[100 100 900 620]);
    plot(rated_stiffness.r_from_blade_root_m, ...
        rated_stiffness.EI_principal_x_Nm2/1e10, ...
        'k-','LineWidth',1.4);
    hold on;
    plot(rated_stiffness.r_from_blade_root_m, ...
        rated_stiffness.EI_principal_y_Nm2/1e10, ...
        'k--','LineWidth',1.4);
    plot(rated_stiffness.r_from_blade_root_m, ...
        rated_stiffness.EI_flap_effective_Nm2/1e10, ...
        'r-','LineWidth',1.4);
    plot(rated_stiffness.r_from_blade_root_m, ...
        rated_stiffness.EI_edge_effective_Nm2/1e10, ...
        'r--','LineWidth',1.4);
    grid on;
    xlabel('Blade span position from root (m)');
    ylabel('Bending stiffness (10^{10} N m^2)');
    title('Baseline and effective bending stiffness');
    legend('Baseline principal x','Baseline principal y', ...
        'Effective flapwise','Effective edgewise','Location','northeast');
    files = [files; export_figure(f,cfg, ...
        '05_effective_bending_stiffness')];
    close(f);
end
end

function files = export_figure(fig,cfg,stem)
files = strings(0,1);
for k = 1:numel(cfg.output.figure_formats)
    format = lower(string(cfg.output.figure_formats{k}));
    filename = fullfile(cfg.paths.figure_dir,stem+"."+format);
    switch format
        case "fig"
            savefig(fig,filename);
        case {"png","pdf"}
            exportgraphics(fig,filename,'Resolution',200);
        otherwise
            error('DTU10MW:FigureFormat','Unsupported figure format %s.',format);
    end
    files(end+1,1) = filename; %#ok<AGROW>
end
end
