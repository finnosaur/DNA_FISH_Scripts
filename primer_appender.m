%% Primer appender
% -depends on MATLAB informatics toolbox
% -requires user to specify desired seqs to append inside of deconvotluion
% primers (ie RT site and T7 promoter)

probe_path = "G:\Finn\20240130_HLB_probe_order\2024-02-07-PaintSHOP-full-probe-file.xlsx";

% table containing your probe set names and sequences
my_probes = readtable(probe_path);

proj_dir = fileparts(probe_path);

% table containing the fwd and rev primers that will be used to deconvolve
% each set (NOTE these are primer sequences, so the rev primer will be rev
% complemented before being appended to your probe sequence)
ps = readtable("G:\Finn\20240130_HLB_probe_order\primer_optimization_v2\subramanian12primerPairs.tsv", 'FileType','text');

fivepr_rt_seq = 'CGTGGTCGCGTCTCA'; 

threepr_t7_seq = 'CCCTATAGTGAGTCGTATTA';


%% get unique targets
targets = unique(my_probes.target);

% boetigger appraoch
% ps_fwd_primers = string(ps_fwd.seq);
% 
% ps_rev_primers = string(ps_rev.seq);

ps_fwd_primers = string(ps.seqSfwd);

ps_rev_primers = string(ps.seqSrev);

%% plot probe locations (must be on same chromosome!)

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
title('HLB Probeset Alignment')
xlabel('chr 6 (bp)')
ylabel('logical')
legend(strcat(string(targets), '--(n=', string(set_size), ')'), "Interpreter","none", "Location","bestoutside")


savefig(fullfile(proj_dir, 'HLB_probe_alignment.fig'))
%% loop thru sets and append

collect_tables = {}; 

for i = 1:numel(targets)
    
    mask = string(my_probes.target) == targets{i};

    disp(targets{i})
    
    % subset our table for current target
    cur_set_t = my_probes(mask, :);

    % get unique outer primers from the paintSHOP list
    cur_fwd_outer_primer = ps_fwd_primers(i);
    cur_rev_outer_primer = ps_rev_primers(i);

    cur_5pr_append_seq = cur_fwd_outer_primer;
    cur_3pr_append_seq = seqrcomplement(cur_rev_outer_primer); % the 3' seq needs to be rev comped

    cur_set_seqs = string(cur_set_t.sequence);

    appended_cur_set_seqs = strcat(cur_5pr_append_seq, fivepr_rt_seq,...
                                cur_set_seqs, threepr_t7_seq, cur_3pr_append_seq);

      % add new column to cur_set_t
    cur_set_t.appended_seqs = appended_cur_set_seqs;
    cur_set_t.app_5pr_inner = repmat(fivepr_rt_seq, size(appended_cur_set_seqs)) ;
    cur_set_t.app_3pr_inner = repmat(threepr_t7_seq, size(appended_cur_set_seqs)) ;



    % % add new column to cur_set_t
    % cur_set_t.appended_seqs = appended_cur_set_seqs;
    cur_set_t.app_5pr = repmat(cur_5pr_append_seq, size(appended_cur_set_seqs)) ;
    cur_set_t.app_3pr = repmat(cur_3pr_append_seq, size(appended_cur_set_seqs)) ;

    % store each new table

    collect_tables{i} = cur_set_t;
    
    disp(i)
    disp('Probe set:' + string(targets{i}))
    disp( strcat('Appending seqs to 5pr of homology region: 5-', cur_5pr_append_seq, '-3',  '5-',fivepr_rt_seq , '-3'))
    disp(strcat('Appending seqs 3pr of homology region: 5-', threepr_t7_seq ,  '-3' , '5-' , cur_3pr_append_seq , '-3'))
    disp('~~~~~~~~~~~~~~')

end

% combine all tables into one (if more than one target)
% initialize table with first probe set
output_table = collect_tables{1}; 

for j = 2:numel(collect_tables)


    output_table = vertcat(output_table, collect_tables{j});

end

%% save




save_name = fullfile(proj_dir, strcat( string(datetime('today')), 'full_probe_table_appended.csv') );

writetable(output_table, save_name, "Delimiter", ",")

