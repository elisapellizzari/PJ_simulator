% debug plots for simulated data

clear all
close all
clc

addpath(genpath(fullfile('..','agata')))

%%

base_folder = 'results';

list_sim = dir(base_folder);
names = {list_sim.name};
valid_idx = contains(names, 'PJ',  'IgnoreCase', true);

list_sim = list_sim(valid_idx);
names = {list_sim.name};
pat = "98";                                                                 % patient for debug

%%

baseline_official_folder = fullfile(base_folder, "DatasetPJ_baseline");

baseline_pats = dir(baseline_official_folder);
b_names = {baseline_pats.name};
b_valid_idx = contains(b_names, 'patient',  'IgnoreCase', true);
b_pats = baseline_pats(b_valid_idx);
b_names = {b_pats.name};

tirs = [];
tars = [];
tbrs = [];
gris = [];
means_g = [];
stds_g = [];
cvs_g = [];
hypos = 0;
for p = 1:length(b_pats)
    data = readtimetable(fullfile(baseline_official_folder, b_names(p)));
    data.Properties.DimensionNames{1} = 'Time';

    tirs = [tirs timeInTarget(data(1:264, :))];
    tars = [tars timeInHyperglycemia(data(1:264, :))];
    tbrs = [tbrs timeInHypoglycemia(data(1:264, :))];
    gris = [gris gri(data(1:264, :))];
    means_g = [means_g meanGlucose(data(1:264, :))];
    stds_g = [stds_g stdGlucose(data(1:264, :))];
    cvs_g = [cvs_g cvGlucose(data(1:264, :))];
    tmp = findHypoglycemicEvents(data(1:264, :));
    tmp = tmp.timeStart;
    if ~isempty(tmp)
        hypos = hypos + 1;
    end

end

disp(prctile(tirs, [25, 50, 75]))
disp(prctile(tars, [25, 50, 75]))
disp(prctile(tbrs, [25, 50, 75]))
disp(hypos/99)
mean(tbrs)
disp(prctile(gris, [25, 50, 75]))
disp(prctile(means_g, [25, 50, 75]))
disp(prctile(stds_g, [25, 50, 75]))
disp(prctile(cvs_g, [25, 50, 75]))


%% plot simuated data

idx_baseline = contains(names, 'baseline',  'IgnoreCase', true);
baseline_sim = list_sim(idx_baseline);
baseline_name = {baseline_sim.name};

dataset_analysis = [];

for p = 1:99
    fig = figure;
    fig.Position = [2700 135 1555 1110];
    for i = 1:length(baseline_name)
        this_sim = baseline_name{i};
        this_pat = fullfile(base_folder, this_sim, "patient_" + p);
        data = readtimetable(this_pat);



        subplot(5, 1, [1:3])
        hold on, grid on
        h1 = plot(data.t, data.glucose, '.-k', LineWidth=2);

        subplot(5, 1, 4)
        hold on
        h2 = stem(data.t, data.cho, "filled", color=[0.3010, 0.7450, 0.9330]);

        subplot(5, 1, 5)
        hold on
        yyaxis left
        h3 = stem(data.t, data.bolus, "filled", color = [0.4660, 0.6740, 0.1880]);
        yyaxis right
        h4 = plot(data.t, data.basal, 'g--');

        axes_list = findall(gcf, 'type', 'axes');
        linkaxes(axes_list, 'x');

    end

    subplot(5,1,[1:3])
    title(sprintf('Glucose - patient %s', num2str(p)))
    ylabel('(mg/dL)')
    yline(180, 'k--')
    yline(70, 'k--')
    legend([h1], {'cgm'})

    subplot(514)
    title('CHO')
    ylabel('(g/min)')
    % --- Annotate CHO labels ---
    labeled_idx = ~ismissing(data.cho_label);
    labeled_ts  = data.t(labeled_idx);
    labeled_vals = data.cho_label(labeled_idx);

    for k = 1:height(data)
        if data.cho(k) > 0
            % Find the nearest labeled timestamp
            [~, nearest_idx] = min(abs(labeled_ts - data.t(k)));
            label = labeled_vals(nearest_idx);

            % Add label above the stem
            text(data.t(k), data.cho(k) + 1, string(label), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',12, ...
                'Color','k');
        end
    end

    subplot(515)
    legend([h3 h4], {'bolus', 'basal'})
    title('Insulin')
    ylabel('Basal (U/min)')
    yyaxis left
    ylabel('Bolus (U/min)')
    xlabel('Time')

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
    % fig.WindowState = 'maximized';

    if ~exist('debug_plots', 'dir')
        mkdir('debug_plots');
    end
    savefig(sprintf('debug_plots/baseline_%s.fig', num2str(p)))
    saveas(fig, sprintf('debug_plots/baseline_%s.png', num2str(p)))

end
close all



%% plot modulations

idx_basal = contains(names, 'basal',  'IgnoreCase', true);
basal_sim = list_sim(idx_basal);
basal_name = {basal_sim.name};

idx_bolus = contains(names, '_bolus',  'IgnoreCase', true);
bolus_sim = list_sim(idx_bolus);
bolus_name = {bolus_sim.name};

idx_meal = contains(names, 'meal',  'IgnoreCase', true);
meal_sim = list_sim(idx_meal);
meal_name = {meal_sim.name};

sim_name_lists = {basal_name, bolus_name, meal_name};
save_names = {'basal', 'bolus', 'meal'};

for s = 1:3
    this_sim_name = sim_name_lists{s};

    fig = figure;
    for i = 1:length(this_sim_name)
        this_sim = this_sim_name{i};
        this_pat = fullfile(base_folder, this_sim, "patient_" + pat);
        data = readtimetable(this_pat);

        if contains(this_sim, 'baseline',  'IgnoreCase', true)
            linewidth=2;
        else
            linewidth=1;
        end

        subplot(6, 1, [1:3])
        hold on, grid on
        plot(data.t, data.glucose, linewidth=linewidth)

        subplot(6, 1, 4)
        hold on
        stem(data.t, data.cho)

        subplot(6, 1, 5)
        hold on
        stem(data.t, data.bolus)

        subplot(6, 1, 6)
        hold on
        plot(data.t, data.basal, linewidth=linewidth)

        axes_list = findall(gcf, 'type', 'axes');
        linkaxes(axes_list, 'x');

    end

    subplot(6,1,[1:3])
    title(sprintf('Glucose - patient %s', pat))
    ylabel('(mg/dL)')
    legend(this_sim_name)
    yline(180, 'k--')
    yline(70, 'k--')

    subplot(614)
    title('CHO')
    ylabel('(g/min)')
    % --- Annotate CHO labels ---
    labeled_idx = ~ismissing(data.cho_label);
    labeled_ts  = data.t(labeled_idx);
    labeled_vals = data.cho_label(labeled_idx);

    for k = 1:height(data)
        if data.cho(k) > 0
            % Find the nearest labeled timestamp
            [~, nearest_idx] = min(abs(labeled_ts - data.t(k)));
            label = labeled_vals(nearest_idx);

            % Add label above the stem
            text(data.t(k), data.cho(k) + 1, string(label), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom', ...
                'FontSize',12, ...
                'Color','k');
        end
    end

    subplot(615)
    title('Bolus')
    ylabel('(U/min)')

    subplot(616)
    title('Basal')
    ylabel('(U/min)')
    xlabel('Time')

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
    fig.WindowState = 'maximized';

    if ~exist('debug_plots', 'dir')
        mkdir('debug_plots');
    end
    savefig(sprintf('debug_plots/%s_modulation_%s.fig', save_names{s}, pat))
    saveas(fig, sprintf('debug_plots/%s_modulation_%s.png', save_names{s}, pat))
    close(fig)

end


%% plot cho/bolus addiction

idx_ht = contains(names, 'hypotreatment',  'IgnoreCase', true);
ht_sim = list_sim(idx_ht);
ht_name = {ht_sim.name};

idx_cbolus = contains(names, 'correctionbolus',  'IgnoreCase', true);
cbolus_sim = list_sim(idx_cbolus);
cbolus_name = {cbolus_sim.name};

sim_name_lists = {ht_name, cbolus_name};
save_names = {'add_ht', 'add_bolus'};

for s = 1:2
    this_sim_name = sim_name_lists{s};

    fig = figure;
    for i = 1:length(this_sim_name)
        this_sim = this_sim_name{i};
        this_pat = fullfile(base_folder, this_sim, "patient_" + pat);
        data = readtimetable(this_pat);

        if contains(this_sim, 'baseline',  'IgnoreCase', true)
            linewidth=2;
        else
            linewidth=1;
        end

        subplot(6, 1, [1:3])
        hold on, grid on
        plot(data.t, data.glucose, linewidth=linewidth)

        subplot(6, 1, 4)
        hold on
        stem(data.t, data.cho)

        subplot(6, 1, 5)
        hold on
        stem(data.t, data.bolus)

        subplot(6, 1, 6)
        hold on
        plot(data.t, data.basal, linewidth=linewidth)

        axes_list = findall(gcf, 'type', 'axes');
        linkaxes(axes_list, 'x');

    end

    subplot(6,1,[1:3])
    title(sprintf('Glucose - patient %s', pat))
    ylabel('(mg/dL)')
    legend(this_sim_name)
    yline(180, 'k--')
    yline(70, 'k--')

    subplot(614)
    title('CHO')
    ylabel('(g/min)')
    % --- Annotate CHO labels ---
    labeled_idx = ~ismissing(data.cho_label);
    labeled_ts  = data.t(labeled_idx);
    labeled_vals = data.cho_label(labeled_idx);

    for k = 1:height(data)
        if data.cho(k) > 0
            % Find the nearest labeled timestamp
            [~, nearest_idx] = min(abs(labeled_ts - data.t(k)));
            label = labeled_vals(nearest_idx);

            % Add label above the stem
            text(data.t(k), data.cho(k) + 1, string(label), ...
                'HorizontalAlignment','center', ...isi
                'VerticalAlignment','bottom', ...
                'FontSize',12, ...
                'Color','k');
        end
    end

    subplot(615)
    title('Bolus')
    ylabel('(U/min)')

    subplot(616)
    title('Basal')
    ylabel('(U/min)')
    xlabel('Time')

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
    fig.WindowState = 'maximized';

    if ~exist('debug_plots', 'dir')
        mkdir('debug_plots');
    end
    savefig(sprintf('debug_plots/%s_%s.fig', save_names{s}, pat))
    saveas(fig, sprintf('debug_plots/%s_%s.png', save_names{s}, pat))
    close(fig)

end