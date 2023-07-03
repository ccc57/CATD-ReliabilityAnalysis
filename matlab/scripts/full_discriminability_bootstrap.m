function discrim_bootstrap = full_discriminability_bootstrap(mats, metric, n_iterations)
% takes two sets of connectomes that are matched and performs
% discriminability
% metric is either 'euclidean' or 'correlation'
% MAYBE ADD INDIVIDUAL MEASURE

% INPUTS %
% mats -- nodeXnodeXsubjXn_measures matrix
% metric -- 'correlation' or 'euclidean'

%%
    n_subj = size(mats,3);
    n_nodes = size(mats,1);
    n_edges = n_nodes*(n_nodes-1)/2;
    n_measures = size(mats,4);

    discrim_bootstrap = 1:n_iterations;
    n_subj_to_randomize = round(n_subj*.8);

    for it = 1:n_iterations
        suborder = randperm(n_subj);
        suborder = suborder(1:n_subj_to_randomize); % here is where I am pulling out the subjs to randomize each time.
        mats_bs = mats(:,:,suborder, :);

        %% convert connectomes to vectors
        vecs = zeros(n_edges,n_measures*n_subj_to_randomize);
        for subj = 1:n_subj_to_randomize
            for meas = 1:n_measures
                curr_mat = mats_bs(:,:,subj,meas);
                curr_vec = curr_mat(find(triu(ones(n_nodes),1)));
                vecs(:,(subj-1)*n_measures+meas) = curr_vec;
            end
        end
        
        %% create subjectwise distance matrix
        dist_mat = pdist2(vecs',vecs',metric);
        
        %% iterate through rows and add to count
        count = 0;
        for subj = 1:n_subj_to_randomize
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
                end
            end
        end
    
        discrim_bootstrap(it) = count/(n_subj_to_randomize*n_measures*(n_subj_to_randomize*n_measures-n_measures)*(n_measures-1));
    end
end