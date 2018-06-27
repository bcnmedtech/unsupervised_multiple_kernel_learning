%% Launch U-MKL 
 
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


%%% LOAD DATA

load Synthetic_velocities.mat

% FEATURES -- cell array of matrices. Each matrix: nvariables x nsamples. 
% all matrices must have the same number of rows

n_features = numel(FEATURES); %%% Number of features
N = size(FEATURES{1},2); %%% Number of samples

% Used to compute the kernel bandwidth, which is calculated feature-wise as the average
% of the pairwise Euclidean distances between each sample and its k-th nearest neighbour
KNN1 = round(N/2); %%% Other commonly used value --> KNN1 = round(sqrt(N));

%%% Number of neighbours used to define the global affinity matrix
KNN2 = round(N/2); %%% Other commonly used value --> KNN2 = round(sqrt(N));

% Further information about the previous parameters can be found in Section 3.2 of
% Sanchez-Martinez_MedIA_2017 (https://doi.org/10.1016/j.media.2016.06.007)
 
% Vector with the feature kind for posterior assignment of the kernel
% type and parameters
    % kind = 1 --> Pattern
    % kind = 2 --> Continuous variable
    % kind = 3 --> Binary variable
    % kind = 4 --> Categorical (ordinal) variable
kind = ones(1,n_features);

for i = 1:n_features
    switch kind(i)
        case 1
            options.Kernel{i}.KernelType = 'exp_l2';
            options.Kernel{i}.Parameters = KNN1;
        case 2
            options.Kernel{i}.KernelType = 'exp_l2';
            options.Kernel{i}.Parameters = KNN1;
        case 3
            options.Kernel{i}.KernelType = 'binary';
            options.Kernel{i}.Parameters = N;
        case 4
            options.Kernel{i}.KernelType = 'ordinal';
            options.Kernel{i}.Parameters = N;
    end
end
options.AffinityNN = KNN2;
       
%%% Unsupervised MKL
[F_data,~,~] = MKL(FEATURES,options);

% Dimensionality reduced space provided by MKL
figure('name','Dimensionality reduced space provided by MKL')
hold on
for iCluster = 1:5
    clustIdx = label==iCluster;
    plot3(F_data(clustIdx,1),F_data(clustIdx,2),F_data(clustIdx,3),'.','MarkerSize',30,...
       'DisplayName',sprintf('Group %i',iCluster));
end

legend('show');
grid on;
xlabel('Dimension 1'); ylabel('Dimension 2'); zlabel('Dimension 3'); 
title('Output space');
hold off; 
   
