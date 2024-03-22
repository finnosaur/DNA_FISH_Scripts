%% DNA FISH PROBE ALIGNMENT VISUALIZER
% Finn Clark, Lionnet Lab, 3/22/2024

% Takes a paintshop probe table as input (can plot multiple probe sets)
% but all sets for a table must be on same chromosome

% give path to painthops table
probe_path = "E:\2023-01-18-PaintSHOP-full-probe-file.txt";


%% Histogram 

% please give title for the project/experiment
myTitle = 'HLB Probe probe Feiyue';
% give chromosome
myChr = '2L';

figure
set_size = zeros(size(targets));
for i = 1:numel(targets)

    cur_target = targets{i};

    mask = string(my_probes.target) == cur_target;

    disp(targets{i})
    
    % subset our table for current target
    cur_set_t = my_probes(mask, :);
    
    histogram(cur_set_t.start, 20)
    % scatter(cur_set_t.start, ones(size(cur_set_t.start)), "filled")
    alpha(0.5)
    hold on

    set_size(i) = numel(cur_set_t.start);

end

title(strcat(myTitle ,'-hist'))
xlabel(strcat(myChr, ' (bp)'))
ylabel('n')
legend(strcat(string(targets), '--(n=', string(set_size), ')'), "Interpreter","none", "Location","bestoutside")

savefig(fullfile(proj_dir, strcat(myTitle, '_hist.fig')))

%% Scatter
figure
set_size = zeros(size(targets));
for i = 1:numel(targets)

    cur_target = targets{i};

    mask = string(my_probes.target) == cur_target;

    disp(targets{i})
    
    % subset our table for current target
    cur_set_t = my_probes(mask, :);

    scatter(cur_set_t.start, ones(size(cur_set_t.start)), "filled")
    alpha(0.5)
    hold on

    set_size(i) = numel(cur_set_t.start);

end
title(myTitle)
xlabel(strcat(myChr, ' (bp)'))
ylabel('logical')
legend(strcat(string(targets), '--(n=', string(set_size), ')'), "Interpreter","none", "Location","bestoutside")


savefig(fullfile(proj_dir, strcat(myTitle, '_scatter.fig')))
disp('saved figures to')
disp(proj_dir)