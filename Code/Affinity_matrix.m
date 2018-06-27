function W = Affinity_matrix( Data , K_NN2 )

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

N = size(Data,2);

[~,IX1] = sort(Data,1,'descend'); %%% columns unchanged
[~,IX2] = sort(Data,2,'descend'); %%% rows unchanged

Wi = ones(N,N);
Wj = ones(N,N);
for i=1:N
    Wi(i,IX1(K_NN2+2:end,i)) = 0;
end
for j=1:N
    Wj(IX2(j,K_NN2+2:end),j) = 0;
end

W = max(Wi,Wj) .* Data;        

end

