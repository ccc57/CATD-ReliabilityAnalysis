% Copyright 2015 Xilin Shen and Corey Horien

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% The individual functional connectome is unique and stable over months to years
% Corey Horien, Xilin Shen, Dustin Scheinost, R. Todd Constable, bioRxiv
% doi: https://doi.org/10.1101/238113 


%this code allows the calculation of which edges are unhelpful in ID (the
%group consistency metric). To calcuate this metric we multiply an edge
%value from time 1 and time 2 across all subjects and compute the
%mean/median for each edge. Edges with high values in this measure are high
%across all individuals in a group and are not helpful in ID.

% load connectivity matrices from all subjects (in the code below,
% "all_se1_org" and "all_se2_orig" are the original matrices).


% all_se1 is obtained from session 1 
% all_se1 is M by M by N matrix, M is the number of nodes, N is the number of subjects 
% all_se2 is M by M by N matrix, from session 2


%%
run_name = 'all';

all_se1_orig = readmatrix(['./ID_scripts/edge_code/data_v21_1/ses1_',run_name,'.csv']).';
all_se2_orig = readmatrix(['./ID_scripts/edge_code/data_v21_1/ses2_',run_name,'.csv']).';

all_se1_org = ones(122, 122, size(all_se1_orig, 2));
all_se2_org = ones(122, 122, size(all_se1_orig, 2));
for i = 1:size(all_se1_orig, 2)

    a = triu(ones(122),1);
    b = triu(ones(122),1);
    a(a > 0) = all_se1_orig(:,i);
    b(b > 0) = all_se2_orig(:,i);
    out1 = (a + a')./(eye(122)+1);
    out2 = (b + b')./(eye(122)+1);
    all_se1_org(:,:,i) = out1;
    all_se2_org(:,:,i) = out2;
end

all_default_se1 = all_se1_org;
all_default_se2 = all_se2_org;
dim = size(all_default_se1);
no_sub = dim(3);
no_node = dim(1);

aa = ones(no_node, no_node);
aa_upp = triu(aa, 1);
upp_id = find(aa_upp);
upp_len = length(upp_id);


all_default_se1 = reshape(all_default_se1, no_node*no_node, dim(3));
all_default_se1 = all_default_se1(upp_id,:);
all_default_se2 = reshape(all_default_se2, no_node*no_node, dim(3));
all_default_se2 = all_default_se2(upp_id,:);

nn_default_s1 = zeros(upp_len, no_sub); 
nn_default_s2 = zeros(upp_len, no_sub); 

% normalize the connectivity matrices for each subject 

for i=1:no_sub;
    c1 = all_default_se1(:,i);
    c1 = c1-mean(c1);
    nn_default_s1(:,i) = c1/norm(c1); 
    
    c2 = all_default_se2(:,i);
    c2 = c2-mean(c2);
    nn_default_s2(:,i) = c2/norm(c2);
end



edge_corr = zeros(no_sub, upp_len); % save the consistency statistics

for i=1:no_sub;
    cur_edge_corr = nn_default_s1(:,i).*nn_default_s2(:,i);
    edge_corr(i,:) = cur_edge_corr;
end

mean_edge_corr = mean(edge_corr);  % mean across subjects
median_edge_corr = median(edge_corr); % median across subjects

sorted_mean_edge_corr = sort(mean_edge_corr,'descend');
sorted_med_edge_corr = sort(median_edge_corr,'descend');


mean_group_con_edge_stat_matrix = zeros(no_node, no_node);
mean_group_con_edge_stat_matrix(upp_id) = mean_edge_corr ;
mean_group_con_edge_stat_matrix = mean_group_con_edge_stat_matrix + mean_group_con_edge_stat_matrix';

med_group_con_edge_stat_matrix = zeros(no_node, no_node);
med_group_con_edge_stat_matrix(upp_id) = median_edge_corr ;
med_group_con_edge_stat_matrix = med_group_con_edge_stat_matrix + med_group_con_edge_stat_matrix';



% threshold the differential statistic matrix and create binary mask for visualization
filename = 'mean_gr_con_edges';

th = [0.001 0.005 0.01 0.02 0.05 0.1];
th_len = round(th*upp_len); 


for i=1:length(th);
    mask = (mean_group_con_edge_stat_matrix>=sorted_mean_edge_corr(th_len(i)));
    %dlmwrite(['./outs_v21_1/',filename, '_',run_name,'_',  num2str(th(i)),'.txt'], mask, ' ');
    clear mask;
end

%  gr_upper = triu(mean_group_con_edge_stat_matrix);
%  flat_gr = gr_upper(gr_upper ~= 0)';

%from here, using the "Connectivity viewer" available here 
% (https://bioimagesuiteweb.github.io/webapp/) allows visualization of
% where the group consistency edges are located in the brain.


