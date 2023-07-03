group = 'mdd'; %mdd, hv, or all
atlas = 'notdBasc122'; %notdSchaefersc or notdBasc122
timepoint = 'v2'; %v4 or v2
mode = 'bootstrap'; % bootstrap, subject, edge, or vanilla
metric = 'euclidean';
iterations = 1000;
reload = false;

tic
if reload || ~strcmp(subfilename, ['references/discriminability/', group, '_for_ml_', timepoint, '.csv'])
    disp("loading matrices")
    subfilename = ['references/discriminability/', group, '_for_ml_', timepoint, '.csv'];
    subfile = readtable(subfilename);
    example = readmatrix(['data/processed/connectomes/sub-20900/ses-v1/sub-20900_ses-v1_task-rest_run-1_space-MNI152NLin2009cAsym_desc-', atlas, '_correlation.csv']);
    subfile.sessionbin = (subfile.session == "v1");
    subfile.subnum = repelem(1:size(subfile, 1)/2, 2)';
        mats = zeros(size(example, 1), size(example, 1), size(subfile, 1)/2, 2);
        for i = 1:size(subfile, 1)
            mats(:, :, subfile.subnum(i),int8(subfile{i, "sessionbin"}) + 1) = readmatrix(strcat('data/processed/connectomes/sub-',string(subfile.subject(i)),'/ses-',string(subfile.old_sess{i}),'/sub-',string(subfile.subject(i)),'_ses-',string(subfile.old_sess{i}),'_task-rest_run-',string(subfile.run(i)),'_space-MNI152NLin2009cAsym_desc-',atlas,'_correlation.csv'));
        end
end
disp(size(mats))
disp("matrices loaded")

if atlas == "notdBasc122"
    mats = mats(1:122, 1:122, :, :); %Remove extraneous dlPFC node
end

if mode == "bootstrap"
    discrim = full_discriminability_bootstrap(mats, metric, iterations);
elseif mode == "subject"
    discrim = full_discriminability(mats, metric, true);
elseif mode == "edge"
    discrim = full_discriminability_edge(mats, metric);
elseif mode == "vanilla"
    discrim = full_discriminability(mats, metric, false);
end

writematrix(discrim, ['data/output/discriminability/',group,'/',group,'_', timepoint, '_', atlas, '_', mode, '_discrim.csv'],'Delimiter',',')
toc