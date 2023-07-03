group = 'mdd';
timepoint = 'v2';
atlas = 'schaefersc';

ICC_table = readtable(['./data/output/icc/',group,'/',group,'_icc_boots_',atlas,'.csv']);
ICC_table = ICC_table(ICC_table.type == "ICC2", :);
%%
ICC_vector = ICC_table(:,"ICC");
ICC_coords = ICC_table(:,"measure");

coords = zeros(size(ICC_coords, 1), 2);
for i = 1:size(ICC_coords, 1)
    coords(i,:) = str2double(strsplit(char(ICC_coords{i,:}), '_'));
end
%%

sort_ICCs = [coords ICC_vector{:,1}];
sort_ICCs = sortrows(sort_ICCs, [1 2]);

%%
sch_mat = zeros(452, 452);
for i = 1:size(sort_ICCs, 1)
    sch_mat(sort_ICCs(i, 1)+1, sort_ICCs(i, 2)+1) = sort_ICCs(i, 3);
    sch_mat(sort_ICCs(i, 2)+1, sort_ICCs(i, 1)+1) = sort_ICCs(i, 3);
end

writematrix(sch_mat, ['./data/output/icc/',group,'/',group,'_',atlas,'_',timepoint,'_icc_matrix.csv'], 'Delimiter',',');

%%
atlas = 'basc';
ICC_table_basc = readtable(['./data/output/icc/',group,'/',group,'_icc_boots_',atlas,'.csv']);
ICC_table_basc = ICC_table_basc(ICC_table_basc.type == "ICC2", :);
%%
ICC_vector_basc = ICC_table_basc(:,"ICC");
ICC_coords_basc = ICC_table_basc(:,"measure");

coords_basc = zeros(size(ICC_coords_basc, 1), 2);
for i = 1:size(ICC_coords_basc, 1)
    coords_basc(i,:) = str2double(strsplit(char(ICC_coords_basc{i,:}), '_'));
end
%%

sort_ICCs_basc = [coords_basc ICC_vector_basc{:,1}];
sort_ICCs_basc = sortrows(sort_ICCs_basc, [1 2]);

%%
basc_mat = zeros(122, 122);
for i = 1:size(sort_ICCs_basc, 1)
    basc_mat(sort_ICCs_basc(i, 1)+1, sort_ICCs_basc(i, 2)+1) = sort_ICCs_basc(i, 3);
    basc_mat(sort_ICCs_basc(i, 2)+1, sort_ICCs_basc(i, 1)+1) = sort_ICCs_basc(i, 3);
end

writematrix(basc_mat, ['./data/output/icc/',group,'/',group,'_',atlas,'_',timepoint,'_icc_matrix.csv'], 'Delimiter',',');