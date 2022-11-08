% Copyright 2015 Xilin Shen and Corey Horien

% This code is released under the terms of the GNU GPL v2. This code
% is not FDA approved for clinical use; it is provided
% freely for research purposes. If using this in a publication
% please reference this properly as: 

% The individual functional connectome is unique and stable over months to years
% Corey Horien, Xilin Shen, Dustin Scheinost, R. Todd Constable, bioRxiv
% doi: https://doi.org/10.1101/238113 


% This code provides a framework for generating confidence intervals
% for functional connectivity-based identification of individual subjects
% via bootstrapping. It assumes an input of pre-calculated correlation
% matrices of dimension NxN for each subject for each session, where N =
% number of nodes in the chosen brain atlas. Each element (i,j) in these
% matrices represents the correlation between the BOLD timecourses of nodes
% i and j during a single fMRI session. In the bootstrapping procedure, a
% subset of the participants are selected each time and ID is performed to
% generate a distribution of ID values. Confidence intervals can then be
% calculated from this distribution.



% Note that this is a companion code to the script
% "ID_example.m"


clear;
clc;


% load connectivity matrices from all subjects (these are named "all_se1_orig"
% and "all_se2_orig" below --> these connectivity matrices must first be vectorized to length M (assuming
% matrices are symmetric this can be done using triu or tril and then
% reshape)




% all_default_se1 is obtained from session 1
% all_default_se1 is M by N matrix, M is the number of edges in the whole connectivity matrix, N is the number of subjects
% all_default_se2 is M by N matrix, from session 2




all_default_se1 = all_se1_orig;
all_default_se2 = all_se2_orig;


N_iteration = 1000; 
rate = zeros(N_iteration,2);

no_sub = size(all_default_se1, 2);

no_sub_to_randomize = round(0.8*no_sub); %pulling out 0.8 of subjs

tic
for it = 1:N_iteration
    
    suborder = randperm(no_sub);
    suborder = suborder(1:no_sub_to_randomize); % here is where I am pulling out the subjs to randomize each time.
    
    all_se1 = all_default_se1(:,suborder);
    all_se2 = all_default_se2(:,suborder);
    
    count1 = 0;
    count2 = 0;
    
    for i=1: no_sub_to_randomize;
        
        % using session 1 as database
        
        tt_corr1 = all_se2(:, i);
        
        tt_to_all1 = corr(tt_corr1, all_se1);
        [~, va_id1] = max(tt_to_all1);
        
        if( i == va_id1)
            count1 = count1+1;
        end
        tt_to_all1_final(i,:) = tt_to_all1;
        % using session 2 as database
        
        tt_corr2 = all_se1(:, i);
        
        tt_to_all2 = corr(tt_corr2, all_se2);
        [~, va_id2] = max(tt_to_all2);
        
        if( i == va_id2)
            count2 = count2+1;
        end  
    end
    rate(it,:) = [count1/no_sub_to_randomize, count2/no_sub_to_randomize];
    
end


 
