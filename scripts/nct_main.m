%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% computes control energy and stability for two state; 
% different choices: can compute control for the whole network and each
% region, or for a subset of regions (e.g. a network such as the rich club)
% default: only global and regional
% Alina Podschun (2024), using code from Urs Braun (2021) and Sebastian Markett (2022)

% required: Urs Braun's network control toolbox: https://github.com/ursbraun/network_control_and_dopamine

% minimum required input: (n-> nodes, s-> participants)
% structural connectomes for each participants (n x n matrix)
% two brain states x0 and xf (n x 1 matrix each for each participant)

% parameter: T (time horizon) default: T=1 
% parameter: rho (penalty term) default: rho=1

% output: 
% cell arrays nct_global, nct_regional and/or nct_subset contain results per participant in each cell
% -> field X contains optimal trajectories
% -> field U contains optimal control energy values
% -> field n contains the error

% matrices nct_global_transf, nct_regional_transf, nct_subset_transf, nct_random_subset_transf contain 
% integrated and transformed values of U, so that columns 1 and 2 represent stability values of state x0 and 
% xf, respectively, and columns 3 and 4 represent control energy needed to transition from x0 -> xf and xf -> x0

% matrix nct_regions_transf_rank contains results of rank-ordering regions according to their NCT metric 
% contribution: The region contributing the most to a metric will have rank 1, the region contributing least will
% have the maximum regional index (equivalent to value of variable num_regions). Rank ordering is done seperately 
% for each participant and NCT metric, and is done on transformed values of U.
% when wanting to correlate regional NCT contribution to any other regional specific on a group level, we strongly 
% advise using this variable for all further analyses.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% adapt: add path to NCT toolbox; define some dataset-specific variables

working_dir = '/path/to/this/NCT-code/';
addpath '/path/to/network_control_and_dopamine-main'

cd(working_dir)

% define number of ROIs in parcellation; we used 219 cortical ROIs
num_regions = 219;

% main data directory, contains subject-specific, ID-named subdirectories
data_dir = '/subpath/to/data';

% output directory
out_dir = '/subpath/to/output_directory';

% load subj IDs
sub_string = fileread('subpath/to/subj_ID_list.txt');
subj_ids = split(sub_string); %List subjects Ids here (folder names)

% connectomes: sub-path to connectome starting from subj-ID folder
conn_path = 'structural/connectome.mat';

% states: sub-paths to start state (x0) and target state (xf)
% e.g. for states in directory /home/data/subj_ID/fuctional/:
x0_path = 'functional/x0.mat';
xf_path = 'functional/xf.mat';

% define NCT parameters:
T = 1; 
rho = 1; 

% define which analysis to do (global & regional, or subset, or both)
% default: only global & regional
glob_reg = 1;
outname_global_regional = 'outname_glob_reg.mat';

% for analysis of role of a subset of regions, set this to 1:
subset = 0;
% and uncomment and adapt following two lines
% outname_subset = 'outname_subset.mat';
% subset_ids = [2 5 7 12 15 19]; % change IDs of regions in subset

% The script will consider these regions a network and calculcate impact on
% control metrics when the network as a whole is prohibited from exerting
% control.

% results of subset analyses can not be directly compared to either global
% or regional analyses because NCT metrics will differ simply because of 
% the nummber of regions prohibited from control. For statistical tests on 
% the role of subnetworks, we recommend comparing results to a null model, 
% e.g. a spin-test based rotation of region IDs. This code contains a .mat 
% of rotated Lausanne atlas IDs (50 rotations). Use the following code to 
% create a null distribution of subset metrics based on those 50 rotations.

rotated_subset = 0;
% for null distribution, set to 1, uncomment and adapt following two lines:
% load([working_dir '/' data_dir '/50_rotated_lausanne250.mat']);
% outname_rotated = 'outname_rotated.mat';

% doesn't have to be based on a spin test if you prefer another mode of
% creating null models; just has to be a num_regions x num_rotations matrix
% named "perm_id" containing rotated region IDs.


%% main part: should need no adaptions
 
% prepare input matrices for and then calculate control metrics

struct = zeros(num_regions,num_regions,length(subj_ids));

for subj=1:length(subj_ids)
    struct(:,:,subj) = importdata([working_dir '/' data_dir '/' subj_ids{subj} '/' conn_path]);
    % stabilize matrix
    A_star(:,:,subj) = struct(:,:,subj)./(eigs(struct(:,:,subj),1)+1)-eye(size(struct(:,:,subj),1));
    A = A_star;
end

x0 = zeros(num_regions,length(subj_ids));
for subj=1:length(subj_ids)
     x0(:,subj) = importdata([working_dir '/' data_dir '/' subj_ids{subj} '/' x0_path]);
end

xf = zeros(num_regions,length(subj_ids));
for subj=1:length(subj_ids)
    xf(:,subj) = importdata([working_dir '/' data_dir '/' subj_ids{subj} '/' xf_path]);
end

if glob_reg == 1
    disp('calculating global and regional metrics')
    [nct_regions nct_global]=nct_analysis_global_regional(A,T,rho,x0,xf);
end

if subset == 1
    disp('calculating subset metrics')
    [nct_subset]=nct_analysis_subset(A,T,rho,x0,xf,subset_ids);
end

if rotated_subset == 1
    disp('calculating nullmodel metrics for subset regions')
    [rotated_subset_members nct_subset_rotated]=nct_analysis_subset_rotated(A,T,rho,x0,xf,subset_ids, perm_id);
end

% integrate and transform values by dividing by t (timesteps)
% needed for statistical comparison, especially between various settings of T

% t is number of timesteps in analysis; transformation calculates mean
% metric per timestep, per subject, and for regional analyses per node
t = (T*1000)+1;

if glob_reg == 1
    % global transformed output variables will have form: subj x metric
    % regional transformed output variables will have form: subj x region x metric
    nct_global_transf = zeros(length(subj_ids),4);
    nct_regions_transf = zeros(length(subj_ids),num_regions,4); 
    nct_regions_transf_rank = zeros(length(subj_ids),num_regions,4); 
    for subj=1:length(subj_ids) 
        nct_global_transf(subj,1:2)=1./(sum(nct_global{subj}.U(:,1:2))/t); 
        nct_global_transf(subj,3:4)=sum(nct_global{subj}.U(:,3:4))/t; 
            for x = 1:219
                nct_regions_transf(subj,x,1) = 1./(sum(nct_regions{subj}.U1(:,x))/t); 
                nct_regions_transf(subj,x,2) = 1./(sum(nct_regions{subj}.U2(:,x))/t); 
                nct_regions_transf(subj,x,3) = sum(nct_regions{subj}.U3(:,x))/t; 
                nct_regions_transf(subj,x,4) = sum(nct_regions{subj}.U4(:,x))/t; 
            end
            % now further transform results; rank order regions seperatedly for
            % participants; nct_regions_transf_rank doesn't contain NCT values, but
            % instead rank of contribution of each region. E.g. if region 2
            % contributed the most to stability of x0 for participant 3,
            % nct_regions_transf_rank(3,2,1) = 1;
            % and if region 10 contributed the least:
            % nct_regions_transf_rank(3,10,1) = 219 
            % 219 in toy data, generally max index of regions
            order_per_region = zeros(num_regions,4);
            a = 'ascend';
            for i = 1:2
                [~,idx] = sort(nct_regions_transf(subj,:,i),a);
                order_per_region(:,i) = idx;
                clear idx
            end
            d = 'descend';
            for i = 3:4
                [~,idx] = sort(nct_regions_transf(subj,:,i),d);
                order_per_region(:,i) = idx;
                clear idx
            end
            rank_per_region = zeros(num_regions,4);
            for i = 1:4
                for r = 1:num_regions
                    rank_per_region(r,i) = find(order_per_region(:,i) == r);
                end
            end
            nct_regions_transf_rank(subj,:,:) = rank_per_region;
    end % end loop over subjects
    savenct=[working_dir '/' out_dir '/' outname_global_regional];
    try
        save(savenct,'nct_regions', 'nct_regions_transf', 'nct_regions_transf_rank', 'nct_global', 'nct_global_transf');
    catch
        disp('failed to save NCT')
    end
end

if subset == 1
    % subset transformed output variables will have form: subj x metric
    for subj=1:length(subj_ids) 
        nct_subset_transf(subj,1:2)=1./(sum(nct_subset{subj}.U(:,1:2))/t); 
        nct_subset_transf(subj,3:4)=sum(nct_subset{subj}.U(:,3:4))/t; 
    end
    savenct=[working_dir '/' out_dir '/' outname_subset];
    try
        save(savenct,'nct_subset', 'nct_subset_transf');
    catch
        disp('failed to save NCT')
    end
end

if rotated_subset == 1
    % rotated subset transformed output variables will be a 1 x subj cell
    % where nct_rotated_transf{1,subj}.U will have form: rotation x metric
    for subj=1:length(subj_ids) 
        % one value per iteration 
        nct_rotated_transf{subj}.U = zeros(size(perm_id,2),4);
        nct_rotated_transf{subj}.U(:,1) = 1./(sum(nct_subset_rotated{subj}.U1,1)/t);
        nct_rotated_transf{subj}.U(:,2) = 1./(sum(nct_subset_rotated{subj}.U2,1)/t);
        nct_rotated_transf{subj}.U(:,3) = sum(nct_subset_rotated{subj}.U3,1)/t;
        nct_rotated_transf{subj}.U(:,4) = sum(nct_subset_rotated{subj}.U4,1)/t;
    end
    savenct=[working_dir '/' out_dir '/' outname_rotated];
    try
        save(savenct,'nct_subset_rotated', 'nct_rotated_transf','rotated_subset_members');
    catch
        disp('failed to save NCT')
    end
end

disp('Complete!')

