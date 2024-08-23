function [rotated_subset_members nct_subset_rotated]=nct_analysis_subset_rotated(A,T,rho,x0,xf,subset_ids,perm_id)
% function computes control energy and stability for two states within an
% task, when a subset of regions is as a network prohibited from control
% Sebastian Markett (2022), 
% using code from Urs Braun (2021)
% adapted by Alina Podschun (2024)

% required: nct scripts, e.g. addpath /home/marketts/tools/network_control_and_dopamine-main

% input: (n-> nodes, s->participants)
% stabilized structural network A (n*n*s matrix)
% state #1: x0 (n*s)
% state #2: xf (n*s)
% subset_ids: 1 x number_subsetregions vector containing IDs of regions in
% subset network
% perm_id = num_regions x num_permutations vector containing permutet (e.g.
% by a spin test) IDs

% parameter: T (time horizon) default: T=1 
% parameter: rho (penalty term) default: rho=1

%% basic info
% how many participants?
nS = size(A,3);
% how many nodes in network?
nN = size(A,2);

%% compute stability & control energy for each node

parfor i=1:nS % loop over participants
    for j=1:size(perm_id,2) % loop over rotated atlases
  
        % set value in identity matrix for subset members = 0; stability and energy
        % when excluding the subset
        D=eye(nN);    
        null_subset = zeros(1,length(subset_ids));
        for f = 1:length(subset_ids)
            k = perm_id(subset_ids(f),j);
            null_subset(f) = k;
        end
        % build matrix for fake subset members to save info which ones have been
        % excluded in which round
        rotated_subset_members{i}.rotation(:,j) = null_subset;
        for r = 1:length(subset_ids)
	        D(null_subset(r), null_subset(r)) = 0;
        end

         % compute stability for condition A->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), x0(:,i), rho,eye(nN));
            % assign to structure 
            nct_subset_rotated{i}.X1(:,j)=trapz(X_opt.^2); nct_subset_rotated{i}.U1(:,j)=trapz(U_opt.^2); nct_subset_rotated{i}.n1(j)=n_err;
        % compute stability for condition B->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), xf(:,i), rho,eye(nN));
            % assign to structure 
            nct_subset_rotated{i}.X2(:,j)=trapz(X_opt.^2); nct_subset_rotated{i}.U2(:,j)=trapz(U_opt.^2); nct_subset_rotated{i}.n2(j)=n_err;
        % compute X U n for condition A->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), xf(:,i), rho,eye(nN));
            % assign to structure 
            nct_subset_rotated{i}.X3(:,j)=trapz(X_opt.^2); nct_subset_rotated{i}.U3(:,j)=trapz(U_opt.^2); nct_subset_rotated{i}.n3(j)=n_err;
         % compute X U n for condition B->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), x0(:,i), rho,eye(nN));
             % assign to structure  
            nct_subset_rotated{i}.X4(:,j)=trapz(X_opt.^2); nct_subset_rotated{i}.U4(:,j)=trapz(U_opt.^2); nct_subset_rotated{i}.n4(j)=n_err;
       
            disp(['participant ' num2str(i) ' at rotation ' num2str(j)])

    end      
end

end
