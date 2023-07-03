function discrim_edge = full_discriminability_edge_old(mats, metric)
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
    n_edges_orig = n_nodes*(n_nodes-1)/2;
    n_edges = (n_nodes*(n_nodes-1)/2) - 1;
    n_measures = size(mats,4);

    discrim_edge = 1:n_edges_orig;

    for it = discrim_edge
        indx = [1:(it-1) (it+1):size(mats, 1)];
        mats_edge = mats(indx, indx, :, :);
        %% convert connectomes to vectors
        vecs = zeros(n_edges,n_measures*n_subj);
        for subj = 1:n_subj
            for meas = 1:n_measures
                curr_mat = mats_edge(:,:,subj,meas);
                curr_vec = curr_mat(find(triu(ones(n_nodes),1)));
                vecs(:,(subj-1)*n_measures+meas) = curr_vec;
            end
        end
        
        %% create subjectwise distance matrix
        dist_mat = pdist2(vecs',vecs',metric);
        
        %% iterate through rows and add to count
        count = 0;
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
                end
            end
        end
    
        discrim_edge(it) = count/(n_subj*n_measures*(n_subj*n_measures-n_measures)*(n_measures-1));
        discrim_edge = zscore(discrim_edge);
    end
end