function [F_data,betas,A] = MKL(FEATURES, options, resume, resume_iterations)

% ================================================================================
%  This Unsupervised Multiple Kernel Learning (U-MKL) package is (c) BCNMedTech
%  UNIVERSITAT POMPEU FABRA
% 
%  This U-MKL package is free software: you can redistribute it and/or modify
%  it under the terms of the GNU Affero General Public License as
%  published by the Free Software Foundation, either version 3 of the
%  License, or (at your option) any later version.
% 
%  This U-MKL package is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU Affero General Public License for more details.
% 
%  You should have received a copy of the GNU Affero General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%  Author:
%  Sergio Sanchez-Martinez 
%
%  Contributors: 
%  Nicolas Duchateau
%  Gemma Piella
%  Constantine Butakoff
% ================================================================================

%function MKL(FEATURES, options)
%FEATURES -- cell array of matrices. Each matrix: npatients x nvariables
%           all matrices must have the same number of rows
%
%options -- structure:
%options.AffinityNN -- number of neighbors of the affinity matrix, by
%                       default sqrt(number of patients). Default = false
%options.NumberOfIterations -- number of iterations to run
%
%resume -- if True, "dump.mat" will be loaded and the iterations resumed
%           The process will restore all the variables
%           including the FEATURES
%resume_iterations -- how many iterations to add to the previous
%           options.NumberOfIterations. This will be added cumulatively on
%           every resume
%
%options.Kernel -- cell array, one cell per one element of FEATURES. 
%   if empty {} - a default 'exp_l2' kernel is used for everything with 5
%   neighbors
%
%   Each cell is:
%   options.Kernel{i}.KernelType -- type of the kernel to use
%   options.Kernel{i}.Parameters -- parameters for the kernel
%
%   Allowed kernels:
%   .KernelType = 'exp_l2' -- exp{-||xi-xj||^2/2sigma^2}
%   .Parameters = [number_of_neighbors] -- number of neighbors to estimate
%                                           sigma, default sqrt(number of patients)

if ~exist('resume','var')
    resume = false;
end

if ~exist('resume_iterations','var')
    resume_iterations = 0;
end

if ~resume

    n_features = numel(FEATURES);
    N = size(FEATURES{1},2);

    if ~exist('options','var')
        options = [];
    end 
    
    %create default options
    if ~isfield(options, 'AffinityNN')
        options.AffinityNN = floor(sqrt(N)); %rule of thumb guess 
    end
    
    if ~isfield(options, 'NumberOfIterations')
        options.NumberOfIterations = 25;
    end
        
    %create default kernel types
    if ~isfield(options, 'Kernel')
        for i=1:numel(FEATURES)
            options.Kernel{i}.KernelType = 'exp_l2';
            options.Kernel{i}.Parameters = floor(sqrt(N)); %rule of thumb guess            
        end 
    end

    %% MKL ALGORITHM

    %%% To compute, or not, S_W_B and S_W_A using MEX files (faster
    %%% computation)
    computeMEX = 1; 

    %% Calculus of the Inputs

    %%% Compute Kernels
    disp('Kernels calculus...');
    KERNELS = zeros(n_features,N,N);
    K_var = zeros(1,length(n_features));
    parfor c=1:n_features 
        fprintf([num2str(c),' ']);
        tmpF = FEATURES{c}';
        tmpF2 = sort(tmpF,'ascend'); %%% For display purposes only
        tmpK = Kernel_Calculus(tmpF,options.Kernel{c}.KernelType, options.Kernel{c}.Parameters,1); %%% a single function only, remove "Kernel_Calculus"
        tmpK2 = Kernel_Calculus(tmpF2,options.Kernel{c}.KernelType, options.Kernel{c}.Parameters,1);
        KERNELS(c,:,:) = tmpK;
        KERNELS2(c,:,:) = tmpK2;
        K_var(c) = var( tmpK(:) );
    end    

    %%% Sum of Kernels, corrected by variance. This amounts to normalizing
    %%% the features by their variance, thus balancing the contribution to
    %%% the neighborhood information encoded in the global affinity matrix
    %%% (W), eq. 2 in Sanchez-Martinez et al.
    %%% (https://doi.org/10.1016/j.media.2016.06.007)
    disp('Kernels sum...');
    tmpM = min(K_var);
    tmpQ = tmpM ./ K_var;
    K_sum  = zeros(N,N);
    parfor c=1:n_features    
        tmpK = Kernel_Calculus(FEATURES{c}',options.Kernel{c}.KernelType, options.Kernel{c}.Parameters,tmpQ(c));
        K_sum = K_sum + tmpK;
    end
    K_sum = K_sum / n_features;

    %%% Affinity matrix: calculated by means of the sum of the Kernels previously extracted. 
    disp('Affinity matrix...');
    W = Affinity_matrix(K_sum, options.AffinityNN); 
    
    %%% Diagonal matrix
    Diag = sum(W,2);

    clear K_sum K_var tmpK tmpK2 tmpF tmpF2 tmpM tmpQ;

    %% Visualize kernels

    figure('units','normalized','position',[0 0 1 1],'name','Visualization of Kernels and Global Affinity Matrix');
    for c=1:(n_features+1)
        subplot(1,n_features+1,c); hold on;
        if c <= n_features
            title(['Kernel ', num2str(c)]);
            imagesc( squeeze(KERNELS2(c,:,:)) );
        else
            title('Global Affinity Matrix');
            imagesc(W);
        end
        axis equal; axis off; axis ij;       
    end

    disp('Starting optimization...');

    clear FEATURES namesT tmpF2 tmpK2 KERNELS2

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Optimization  

    %%%-----------Optimizacion of A (projection matrix)--------------

    betas = ones(n_features,1)/n_features;      % Kernel weigths
    max_it = options.NumberOfIterations;        % Max. number of iterations in the optimization
    conv_crit = 1e-6;                           % Convergence criterion
    gap = 0;                                    % Minimization value in LowSpace
    gap_old = Inf;                              % Previous value of minimization
    t = 1; 

    ZERO = 10^-6;

    %%% Definition of the variables to save at each iteration
    energy_complete = NaN(max_it,2);  % Value of the final minimization
    min_criteriaA = NaN(max_it,2);    % Value of the first minimization
    min_criteriaB = NaN(max_it,2);    % Value of the second minimization
    min_criteriab = NaN(max_it,2);    % Value of the second minimization
    matrix_rank = NaN(max_it,1);

    K_tot_ALL = permute(KERNELS,[3,2,1]);

    clear KERNELS;
    
    save iteration_preload;
else %if ~resume     
    load iteration_preload;
    load iteration_progress;
    max_it = max_it + resume_iterations;
end   
    

%Start: For the parallelization
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    nworkers = 1;
else
    nworkers = poolobj.NumWorkers*3; % Each worker to run 3 jobs
end


if nworkers~=1
    values = create_ivalues(nworkers, N); % The arrays will start from 0
else 
    values = [0,N];
end

start_values = values(1:end-1);
end_values = values(2:end);

while (t <= max_it) && (abs(gap_old-gap) > conv_crit)
    
    clc;
    disp('*********************');
    fprintf('New CVX iteration: %d / %d \n',t,max_it);
    disp('*********************');
    
    gap_old = gap;
    
    %%%-----------Optimizacion of A (projection matrix)--------------
    
    clear A;
            
    fprintf('computeSWB... [%d/%d]\n',t,max_it);
    
    S_W_B = zeros(N);
    S_D_B = zeros(N);
    
    if computeMEX == 1

        parfor i=1:numel(start_values)
            [S_W_B1 , S_D_B1] = computeSWB(K_tot_ALL, betas, W, Diag,...
                start_values(i), end_values(i) ); %%% MEX-file
            S_W_B = S_W_B + S_W_B1;
            S_D_B = S_D_B + S_D_B1;
        end

        %%% Fill in below diagonal
        S_W_B = S_W_B + triu(S_W_B,1)'; 
        S_D_B = S_D_B + triu(S_D_B,1)';

        %%% To guarantee that the matrices are symmetric
        S_W_B = 0.5 * (S_W_B + S_W_B'); 
        S_D_B = 0.5 * (S_D_B + S_D_B');
        
    else
        
        fprintf('NONMEXcomputeSWB... [%d/%d]\n',t,max_it);
        nonmexS_W_B = zeros(N,N);
        nonmexS_D_B = zeros(N,N);     
        for i=1:N
            Ki = squeeze(K_tot_ALL(i,:,:));
            for j=(i+1):N %%% W is symmetric
                if(W(i,j) ~= 0) %%% W is sparse
                    K_tot = Ki - squeeze(K_tot_ALL(j,:,:));
                    tmp = K_tot*betas;
                    nonmexS_W_B = nonmexS_W_B + 2 * W(i,j)*(tmp*tmp'); %%% sum twice
                end
            end        
            tmp = Ki*betas;
            nonmexS_D_B = nonmexS_D_B + Diag(i)*(tmp*tmp');
        end   
        %%% To guarantee that the matrices are symmetric
        nonmexS_W_B = 0.5 * (nonmexS_W_B + nonmexS_W_B'); 
        nonmexS_D_B = 0.5 * (nonmexS_D_B + nonmexS_D_B');
        S_W_B = nonmexS_W_B;
        S_D_B = nonmexS_D_B;

    end
            
    % Check if S_W_B & S_D_B are Positive Definite Matrices (PDM). If not,
    % the eig function will return complex values. They will be positive
    % definite if all their eigenvalues are greater than zero. The following 
    % lines convert these matrices to PDM
    disp('To_PDM...');
    tmp = eye(N,N) * ZERO;
    S_W_B = S_W_B + tmp;
    S_D_B = S_D_B + tmp;
    
    %%% Generalized eigenvalue decomposition problem
    opts.isreal = 1;
    opts.issym = 1;
    disp('eigs...');
    matrix_rank(t,1) = rank(S_D_B*n_features);
    [V,D] = eigs(S_W_B,S_D_B*n_features,N,'sa',opts); 

    %%% Arrange the eigenvalues and associated eigenvectos (D and V) 
    %%% from smallest to largest
    if D(end,end) == max(max(D))
        D = diag(D);
    else
        V = fliplr(V);
        D = diag(D); D = flipud(D);
    end
        
    %%% Projection matrix
    A = V; 
    A = A(:,D>ZERO); %%% remove components associated to eig = 0
       
    %%% Calculus of the expresion to minimize
    min_criteriaA(t,1) = trace(A'*S_W_B*A);
    min_criteriaA(t,2) = trace(A'*S_D_B*A);
    
    %%%-----------Optimizing Betas (Kernel weights)--------------
    
    % The parallelized computation of SWA is slower than the non-parallelized
    % version when the number of features is high (suggestion, higher than 
    % 50 features), due to data transfer overhead. 
    
    computeMEX = 1;
    
    if computeMEX == 1
    
        fprintf('computeSWA... [%d/%d]\n',t,max_it);

        S_W_A = zeros(n_features);
        S_D_A = zeros(n_features);
        parfor i=1:numel(start_values)   
            [S_W_A1 , S_D_A1] = computeSWA(K_tot_ALL, A, W, Diag,...
                start_values(i), end_values(i) ); %%% MEX-file

            S_W_A = S_W_A + S_W_A1;
            S_D_A = S_D_A + S_D_A1;
        end
        S_W_A = S_W_A + triu(S_W_A,1)'; %%% Fill in below diagonal
        S_D_A = S_D_A + triu(S_D_A,1)';


        S_W_A = 0.5 * (S_W_A + S_W_A'); 
        S_D_A = 0.5 * (S_D_A + S_D_A');
        
    else
    
        fprintf('NONMEXcomputeSWA... [%d/%d]\n',t,max_it);
        nonmexS_W_A = zeros(n_features,n_features);
        nonmexS_D_A = zeros(n_features,n_features);     
        for i=1:N
            Ki = squeeze(K_tot_ALL(i,:,:));
            for j=(i+1):N %%% W is symmetric
                if(W(i,j) ~= 0) %%% W is sparse
                    K_tot = Ki - squeeze(K_tot_ALL(j,:,:));
                    tmp = K_tot'*A;
                    nonmexS_W_A = nonmexS_W_A + 2 * W(i,j)*(tmp*tmp'); %%% sum twice
                end
            end        
            tmp = Ki'*A;
            nonmexS_D_A = nonmexS_D_A + Diag(i)*(tmp*tmp');
        end             

    nonmexS_W_A = 0.5 * (nonmexS_W_A + nonmexS_W_A'); 
    nonmexS_D_A = 0.5 * (nonmexS_D_A + nonmexS_D_A');

    S_W_A = nonmexS_W_A;
    S_D_A = nonmexS_D_A;
    
    end
        
    % CVX toolbox: semidefinite programming to optimize betas (see Lin_PAMI_2011
    % for further information)
    clear betas B;
    
    disp('CVX...');
    
    cvx_begin 
    
        %%% With a high number of features and samples, high precision may be
        %%% needed to reach the same results at different launches
        %cvx_precision high 
        variable B(n_features,n_features) symmetric
        variable betas(n_features,1) 
        minimize ( trace(S_W_A * B) )
        subject to     
            trace( S_D_A * B ) == 1
            for m=1:n_features
                betas(m) >= 0;
            end
            [1 , betas' ; betas , B] == semidefinite(n_features+1);                        
    cvx_end 

    clc;
    disp('*********************');
    fprintf('Ongoing CVX iteration: %d / %d \n',t,max_it);
    disp('*********************');

    constant = sum(betas);
    %%% force sum of betas = 1
    betas = betas / constant; 
    
    % Computation of minimization criteria
    min_criteriaB(t,1) = trace(S_W_A*B); 
    min_criteriaB(t,2) = trace(S_D_A*B);
        
    min_criteriab(t,1) = betas'*S_W_A*betas;
    min_criteriab(t,2) = betas'*S_D_A*betas;

    %%% Cost function evaluation with the refreshed parameters (A & betas)   
    fprintf('computeENERGY... [%d/%d]\n',t,max_it);
    
    computeMEX = 1;
    
    %%% Compute the global energy. Eq.4 in Sanchez-Martinez et al. 
    %%% (https://doi.org/10.1016/j.media.2016.06.007)
    if(computeMEX == 1)

        gap = 0;
        constr = 0;
        parfor i=1:numel(start_values)      
            [gap1,constr1] = computeENERGY(K_tot_ALL,betas,A,W,Diag,...
                start_values(i), end_values(i) );

            gap = gap+gap1;
            constr = constr+constr1;
        end
        
    else
    
        fprintf('NONMEXcomputeENERGY... [%d/%d]\n',t,max_it);
        NONMEXgap = 0;     
        NONMEXconstr = 0;     
        for i=1:N        
            Ki = squeeze(K_tot_ALL(i,:,:));
            for j=(i+1):N %%% W is symmetric
                if(W(i,j) ~= 0) %%% W is sparse
                    Kj = squeeze(K_tot_ALL(j,:,:));
                    NONMEXgap = NONMEXgap + 2 * sum ( ( A'*(Ki-Kj)*betas ).^2 ) * W(i,j); %%% sum twice
                end
            end        
            NONMEXconstr = NONMEXconstr + sum ( ( A'*Ki*betas ).^2 ) * Diag(i);
        end  
        
        gap = NONMEXgap;
        constr = NONMEXconstr;
   
    end
    
    energy_complete(t,1) = gap; 
    energy_complete(t,2) = constr; 
    
    if t == 2
        conv_crit = abs(gap_old-gap)*0.02;
    end
    
    t = t+1;     
    save iteration_progress energy_complete min_criteriab min_criteriaA min_criteriaB betas matrix_rank t max_it gap_old gap;
   
end

%%% Clear variables
clear tmp K_tot Ki Kj W Diag B D V;
if(computeMEX == 1)
    clear S_D_A S_D_B S_W_A S_W_B;
else
    clear nonmexS_D_A nonmexS_D_B nonmexS_W_A nonmexS_W_B;
end

%%% Linear combination of Kernels
KERNELS = permute(K_tot_ALL,[3,2,1]);
clear K_tot_ALL;
Kernel = zeros(N,N);
parfor c=1:n_features
    Kernel = Kernel + squeeze(KERNELS(c,:,:))*betas(c);
end

%%% Data projected in the learned representation space
F_data = zeros(N,size(A,2));    
parfor i =1:N
    F_data(i,:) = (Kernel(i,:)*A)';
end


%% Final Plot

% Choose Output Space dimensions to plot
dim1 = 1;
dim2 = 2;
dim3 = 3;

% Global Energy 
figure('units','normalized','position',[0 0 1 1],'name','MKL Optimization');
subplot(2,3,1);
hold on; grid on;
plot(energy_complete(:,1),'b','Marker','+','MarkerSize',10);
plot(min_criteriab(:,1),'r','Marker','.');
legend('Global energy','Energy betas');
axis square;

% Projection matrix --- Energy and constraint 
subplot(2,3,2);
hold on; grid on;
plot(min_criteriaA(:,1),'r','Marker','+');
plot(min_criteriaA(:,2),'k','Marker','.');
legend('Energy proj. mat.','Constraint proj. mat.');
axis square;

% Contraint B
subplot(2,3,3);
hold on; grid on;
plot(min_criteriaB(:,2),'k','Marker','.');
legend('Constraint B');
axis square;

% Constraints --- Global and betas
subplot(2,3,4);
hold on; grid on;
plot(energy_complete(:,2),'b','Marker','+','MarkerSize',10);
plot(min_criteriab(:,2),'r','Marker','.');
legend('Global constraint','Constraint betas');
axis square;

% Plot weights (betas)
subplot(2,3,5);
hold on; grid on;
bar(betas);
title('Weights'),
axis square;

% Plot Output Space
subplot(2,3,6);
hold on; grid on;
scatter(F_data(:,1),F_data(:,2),100,'k','filled'); hold on;
xlabel('Dim1'); ylabel('Dim2'); zlabel('Dim3');
title('Output space');
axis square;

save dump;

disp('Finished');
