clear all

filename = ['MSiCOR_DATA_after_2006_10_codes_have_covariates.csv'];
MAT = csvread(filename);

filename = ['Cluster_membership_ordered.csv'];
MS_patient_membership = csvread(filename);

N = size(unique(MAT(:,1)),1);

sequences_temp = cell(1,N);

patient_num = 1;
sequences_temp{1} = MAT(1,2);
seq_lengths = zeros(N,1);

for i = 2:size(MAT,1)
    if(MAT(i,1) == MAT(i-1,1))
        sequences_temp{patient_num} = [sequences_temp{patient_num} MAT(i,2)];
    else
        patient_num = patient_num + 1;
        sequences_temp{patient_num} = MAT(i,2);
    end
end

number_of_changes = zeros(N,1);
seq_lengths = zeros(N,1);
rate_of_changes_per_year = zeros(N,1);

for i = 1:N
    seq_lengths(i) = size(sequences_temp{i},2);
    num_shifts = 0;
    for j = 2:seq_lengths(i)
        if(sequences_temp{i}(j) ~= sequences_temp{i}(j-1))
           num_shifts = num_shifts + 1;
        end
    end
    number_of_changes(i) = num_shifts;
    rate_of_changes_per_year(i) = 4*number_of_changes(i)/(seq_lengths(i)-1);
end

clus_1_values = rate_of_changes_per_year(MS_patient_membership(:,2) == 1);
clus_2_values = rate_of_changes_per_year(MS_patient_membership(:,2) == 2);
clus_3_values = rate_of_changes_per_year(MS_patient_membership(:,2) == 3);

length(clus_1_values)
length(clus_2_values)
length(clus_3_values)

nan_positions_clus_1 = find(isnan(clus_1_values));
nan_positions_clus_2 = find(isnan(clus_2_values));
nan_positions_clus_3 = find(isnan(clus_3_values));

clus_1_values_refined = clus_1_values;
clus_1_values_refined(nan_positions_clus_1) = [];

clus_2_values_refined = clus_2_values;
clus_2_values_refined(nan_positions_clus_2) = [];

clus_3_values_refined = clus_3_values;
clus_3_values_refined(nan_positions_clus_3) = [];


clusterwise_rate = [mean(clus_1_values_refined), mean(clus_2_values_refined), mean(clus_3_values_refined)];
clusterwise_rate_se = [std(clus_1_values_refined)/sqrt(length(clus_1_values_refined)),...
                       std(clus_2_values_refined)/sqrt(length(clus_2_values_refined)),...
                       std(clus_3_values_refined)/sqrt(length(clus_3_values_refined))];

filename = ['num_of_DMT_changes_per_year.csv'];
csvwrite(filename,clusterwise_rate)

filename = ['SE_num_of_DMT_changes_per_year.csv'];
csvwrite(filename,clusterwise_rate_se)
