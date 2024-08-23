function [nct_regions nct_global]=nct_analysis_global_regional(A,T,rho,x0,xf)
% function computes control energy and stability for two states within an
% task; fist for the whole network, then after iteratively prohibiting each
% node from exhibiting control
% Sebastian Markett (2022), 
% using code from Urs Braun (2021)
% adapted by Alina Podschun (2024)

% required: nct scripts, e.g. addpath /home/marketts/tools/network_control_and_dopamine-main

% input: (n-> nodes, s->participants)
% stabilized structural network A (n*n*s matrix)
% state #1: x0 (n*s)
% state #2: xf (n*s)

% parameter: T (time horizon) default: T=1 
% parameter: rho (penalty term) default: rho=1

%% basic info
% how many participants?
nS = size(A,3);
% how many nodes in network?
nN = size(A,2);

%% compute stability & control energy for each node

parfor i=1:nS % loop over participants
    for j=1:nN % loop over regions

        D=eye(nN); D(j,j)=0; % iterativly remove one node
        
        % compute stability for condition A->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), x0(:,i), rho,eye(nN));
            % assign to structure 
            nct_regions{i}.X1(:,j)=trapz(X_opt.^2); nct_regions{i}.U1(:,j)=trapz(U_opt.^2); nct_regions{i}.n1(j)=n_err;
        % compute stability for condition B->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), xf(:,i), rho,eye(nN));
            % assign to structure 
            nct_regions{i}.X2(:,j)=trapz(X_opt.^2); nct_regions{i}.U2(:,j)=trapz(U_opt.^2); nct_regions{i}.n2(j)=n_err;
        % compute X U n for condition A->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), xf(:,i), rho,eye(nN));
            % assign to structure 
            nct_regions{i}.X3(:,j)=trapz(X_opt.^2); nct_regions{i}.U3(:,j)=trapz(U_opt.^2); nct_regions{i}.n3(j)=n_err;
         % compute X U n for condition B->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), x0(:,i), rho,eye(nN));
             % assign to structure  
            nct_regions{i}.X4(:,j)=trapz(X_opt.^2); nct_regions{i}.U4(:,j)=trapz(U_opt.^2); nct_regions{i}.n4(j)=n_err;

            disp(['regional analysis: participant ' num2str(i) ' at node ' num2str(j)])

    end      
end
  
%% compute stability & control energy for whole network

parfor i=1:nS % loop over participants

        D=eye(nN); % control input to all nodes
        
          % compute stability for condition A->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), x0(:,i), rho,eye(nN));
         nct_global{i}.X(:,1)=trapz(X_opt.^2); nct_global{i}.U(:,1)=trapz(U_opt.^2); nct_global{i}.n(1)=n_err;
          % compute stability for condition B->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), xf(:,i), rho,eye(nN));
         nct_global{i}.X(:,2)=trapz(X_opt.^2); nct_global{i}.U(:,2)=trapz(U_opt.^2); nct_global{i}.n(2)=n_err;
        % compute energy for condition A->B
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, x0(:,i), xf(:,i), rho,eye(nN));
        % assign to structure 
        nct_global{i}.X(:,3)=trapz(X_opt.^2); nct_global{i}.U(:,3)=trapz(U_opt.^2); nct_global{i}.n(3)=n_err;
         % compute energy for condition B->A
        [X_opt, U_opt, n_err] = optim_fun(A(:,:,i), T, D, xf(:,i), x0(:,i), rho,eye(nN));
         nct_global{i}.X(:,4)=trapz(X_opt.^2); nct_global{i}.U(:,4)=trapz(U_opt.^2); nct_global{i}.n(4)=n_err;
   
end 

end