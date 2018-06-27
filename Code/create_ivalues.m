function ivalues = create_ivalues(n_splits, N)

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

i = N:-1:1;
a = cumsum(i)';
t = ceil(a(end)/n_splits);
ts = t*(1:n_splits);
mask = sum(bsxfun(@le,a,ts),2);
dmask = [-1;diff(mask)];
ivalues = [find(dmask==-1)-1;N]; %the arrays will start from 0 in C++, but last must be N, see the MEX code
