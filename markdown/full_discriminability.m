function discrim = full_discriminability(mats,metric, n_iterations)
% takes two sets of connectomes that are matched and performs
% discriminability
% metric is either 'euclidean' or 'correlation'
% MAYBE ADD INDIVIDUAL MEASURE

% INPUTS %
% mats -- nodeXnodeXsubjXn_measures matrix
% metric -- 'correlation' or 'euclidean'

%%
metric = 'euclidean';
n_subj = size(mats,3);
n_nodes = size(mats,1);
n_edges = n_nodes*(n_nodes-1)/2;
n_measures = size(mats,4);

n_iterations = 10; 
if
discrim = 1:n_iterations;
n_subj_to_randomize = round(n_subj*.8);

for it = 1:n_iterations
    
    suborder = randperm(n_subj);
    suborder = suborder(1:n_subj_to_randomize); % here is where I am pulling out the subjs to randomize each time.

    mats = mats(:,:,suborder, :);

    %% convert connectomes to vectors
    vecs = zeros(n_edges,n_measures*n_subj);
    for subj = 1:n_subj
        for meas = 1:n_measures
            curr_mat = mats(:,:,subj,meas);
            curr_vec = curr_mat(find(triu(ones(n_nodes),1)));
            vecs(:,(subj-1)*n_measures+meas) = curr_vec;
        end
    end
    
    %% create subjectwise distance matrix
    dist_mat = pdist2(vecs',vecs',metric);
    
    %% iterate through rows and add to count
    count = 0;
    subjcounts = zeros(n_subj, 1);
    for subj = 1:n_subj
        dist_inds = (subj-1)*n_measures + (1:n_measures);
        for row = dist_inds
            full_row = dist_mat(row,:);
            cross_scan_inds = dist_inds;
            cross_scan_inds(dist_inds==row) = [];
            for cross_scan_ind = cross_scan_inds
                curr_comparison = full_row(cross_scan_ind);
                other_comparisons = full_row;
                other_comparisons(dist_inds) = [];
                count = count + length(find(other_comparisons>=curr_comparison));
                subjcounts(subj) = subjcounts(subj) + length(find(other_comparisons>=curr_comparison));
            end
        end
    end

    discrim = count/(n_subj*n_measures*(n_subj*n_measures-n_measures)*(n_measures-1));
    subjdiscrim = subjcounts/((n_subj*n_measures-n_measures)*(n_measures-1)*n_measures);
writematrix(subjdiscrim, 'data/output/mdd_subj_discrim.csv','Delimiter',',')

end


subfile = readtable('references/mdd_for_ml.csv');
example = readmatrix('data/processed/connectomes/sub-20900/ses-v1/sub-20900_ses-v1_task-rest_run-1_space-MNI152NLin2009cAsym_desc-notdSchaefersc_correlation.csv');
subfile.sessionbin = (subfile.session == "v1");
subfile.subnum = repelem(1:size(subfile, 1)/2, 2)';
mats = zeros(size(example, 1), size(example, 1), size(subfile, 1)/2, 2);
for i = 1:size(subfile, 1)
    mats(:, :, subfile.subnum(i),int8(subfile{i, "sessionbin"}) + 1) = readmatrix(strcat('data/processed/connectomes/sub-',string(subfile.subject(i)),'/ses-',string(subfile.old_sess{i}),'/sub-',string(subfile.subject(i)),'_ses-',string(subfile.old_sess{i}),'_task-rest_run-',string(subfile.run(i)),'_space-MNI152NLin2009cAsym_desc-notdSchaefersc_correlation.csv'));
end