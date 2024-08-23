function [nct_subset]=nct_analysis_subset(A,T,rho,x0,xf,subset_ids)
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

% parameter: T (time horizon) default: T=1 
% parameter: rho (penalty term) default: rho=1


%% basic info
% how many participants?
nS = size(A,3);
% how many nodes in network?
nN = size(A,2);

% set value in identity matrix for subset members = 0; stability and energy
% when excluding the subset
D=eye(nN);
for r = 1:length(subset_ids)
	D(subset_ids(r), subset_ids(r)) = 0;
end
  
%% compute stability & control energy for whole network without rc

parfor i=1:nS % loop over participants
        
         % compute stability for condition A->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), x0(:,i), rho,eye(nN));
         nct_subset{i}.X(:,1)=trapz(X_opt.^2); nct_subset{i}.U(:,1)=trapz(U_opt.^2); nct_subset{i}.n(1)=n_err;
          % compute stability for condition B->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), xf(:,i), rho,eye(nN));
         nct_subset{i}.X(:,2)=trapz(X_opt.^2); nct_subset{i}.U(:,2)=trapz(U_opt.^2); nct_subset{i}.n(2)=n_err;
        % compute energy for condition A->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), xf(:,i), rho,eye(nN));
        % assign to structure 
        nct_subset{i}.X(:,3)=trapz(X_opt.^2); nct_subset{i}.U(:,3)=trapz(U_opt.^2); nct_subset{i}.n(3)=n_err;
         % compute energy for condition B->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), x0(:,i), rho,eye(nN));
         nct_subset{i}.X(:,4)=trapz(X_opt.^2); nct_subset{i}.U(:,4)=trapz(U_opt.^2); nct_subset{i}.n(4)=n_err;
  
end
  

end
