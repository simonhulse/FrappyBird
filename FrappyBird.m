% FrappyBird.m
% Simon Hulse
% simonhulse@protonmail.com
% Last Edited: Wed 26 Feb 2025 05:59:11 PM EST

function FrappyBird
    disp(' ____  ____    __    ____  ____  _  _    ____  ____  ____  ____  ');
    disp('( ___)(  _ \  /__\  (  _ \(  _ \( \/ )  (  _ \(_  _)(  _ \(  _ \ ');
    disp(' )__)  )   / /(__)\  )___/ )___/ \  /    ) _ < _)(_  )   / )(_) )');
    fprintf('(__)  (_)\\_)(__)(__)(__)  (__)   (__)   (____/(____)(_)\\_)(____/ \n\n');

    [~, parent, ~] = fileparts(pwd) ;
    fprintf('Processing data in the directory %s\n\n', parent);

    listing = dir('*.csv');
    n_files = length(listing);
    if n_files == 0
        error('No CSV files found in the directory!');
    end

    root_dir = 'output';
    [~, ~, ~] = mkdir(root_dir);

    fraps = configureDictionary('string', 'cell');
    for fileinfo = listing'
        csvfile = fileinfo.name;

        tokens = regexp(csvfile, '(.*) ([A-Z]).*\.csv', 'tokens');
        conditions = tokens{1}{1};
        letter = tokens{1}{2};

        conditions_dir = sprintf('%s/%s', root_dir, conditions');
        letter_dir = sprintf('%s/%s', conditions_dir, letter');

        fig1_path = sprintf('%s/1.pdf', letter_dir);
        fig2_path = sprintf('%s/2.pdf', letter_dir);
        frap_path = sprintf('%s/frap.mat', letter_dir);
        stats_path = sprintf('%s/stats.txt', letter_dir);

        fprintf('Processing %s\n', csvfile);
        [~, ~, ~] = mkdir(conditions_dir);
        [~, ~, ~] = mkdir(letter_dir);

        % Load csv file
        data = table2array(readtable(csvfile, 'NumHeaderLines', 1));
        index = data(:, 1);
        mean1 = data(:, 3);
        mean2 = data(:, 7);
        mean3 = data(:, 11);

        % Substract background from bleached and non-bleached spots
        roi1 = mean1 - mean3;
        roi2 = mean2 - mean3;
        avgroi1_prb = mean(roi1(1:9), "all");
        avgroi2_prb = mean(roi2(1:9), "all");

        % Convert frames to time points
        time = (7.0 / 40.0) * index;
        % Average roi1 intensity prebleach to normalize
        I_normal1 = roi1 / avgroi1_prb;
        % Average roi2 intensity prebleach to normalize
        I_normal2 = roi2 / avgroi2_prb;

        frap = I_normal1 ./ I_normal2;
        save(frap_path, 'frap');
        fprintf('--> Saved FRAP vector to %s\n', frap_path);

        if isKey(fraps, conditions_dir)
            x = double(char(letter)) - 64;
            old_avg = lookup(fraps, conditions_dir);
            new_avg = ((x - 1) * old_avg{1} + frap) / x;
            fraps = insert(fraps, conditions_dir, {new_avg}, Overwrite=true);
        else
            fraps = insert(fraps, conditions_dir, {frap});
        end

        time_slice = time(11:end);
        frap_slice = frap(11:end);
        [fit_result, ~] = tau_fit(time_slice, frap_slice);

        I = fit_result.I; a = fit_result.a; b = fit_result.b; g = fit_result.g;
        mf = (a + g) / (1.0 - (I - a - g));
        t_half = log(2.0) / b;

        % Determine errors
        values = coeffvalues(fit_result);
        confints = confint(fit_result, 0.6827);
        stdevs = num2cell(values - confints(1, :));
        [I_err, a_err, b_err, d_err, g_err] = stdevs{:};
        mf_error = 0  %% TODO: compute
        t_half_error = (log(2) / b^2) * b_err

        text = sprintf('mf: %6f ± %6f\nt_half: %6f ± %6f', mf, mf_error, t_half, t_half_error);
        fileID = fopen(stats_path, 'w');
        fprintf(fileID, text);
        fclose(fileID);
        fprintf('%s\n', text);
        fprintf('--> Saved `mf` and `nt_half` to %s\n', stats_path);

        % Plot normalized data
        fig1 = figure;
        plot(time, frap, 'r.');
        xlabel('Time (s)');
        ylabel('Norm. Intensity');
        exportgraphics(fig1, fig1_path);
        close(fig1);
        fprintf('--> Saved FRAP figure to %s\n', fig1_path);

        fig2 = figure('Name', 'Double fit');

        % Plot fit
        subplot(2, 1, 1);
        h = plot(fit_result, time_slice, frap_slice);
        legend(h, 'FRAP', 'Fit', 'Location', 'SouthEast');
        ylabel('Norm. Intensity');
        grid on;

        % Plot residuals.
        subplot(2, 1, 2);
        h = plot(fit_result, time_slice, frap_slice, 'residuals');
        xlabel('Time (s)');
        ylabel('Residual');
        grid on;

        exportgraphics(fig2, fig2_path);
        close(fig2);
        fprintf('--> Saved fit figure to %s\n\n', fig2_path);
    end

    fprintf('Saving FRAP averages...\n');
    keys_vals = entries(fraps);
    for i = 1:height(keys_vals)
        directory = keys_vals{i, 1};
        frap_avg = keys_vals{i, 2}{1};
        frap_avg_path = sprintf('%s/frap_avg.mat', directory);
        save(frap_avg_path, 'frap_avg');
        fprintf('--> Saved FRAP average to %s\n', frap_avg_path);
    end

    fprintf('\n    ---------\n');
    disp('   ( DONE!!! )');
    disp('    ---------');
    disp('       /');
    disp('___()>');
    disp('\___)');
end

function [fit_result, gof] = tau_fit(x, y)
    ft = fittype('I-a*exp(-b*x)-g*exp(-d*x)', 'independent', 'x', 'dependent', 'y');

    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.StartPoint = [0.622475086001227 0.587044704531417 0.207742292733028 0.301246330279491 0.470923348517591];

    [fit_result, gof] = fit(x, y, ft, opts);
end
