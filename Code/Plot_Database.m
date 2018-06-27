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

groups = unique(label);
color = 'brgmk';
figure;
% Plot curves
for i = 1:numel(groups)
    subplot(numel(groups),1,i);
    plot(curves(:,label==groups(i)),color(i)); grid on; box on; axis tight;
    hold on;
    %%% vertical bars demarking different segments
    plot([500,500],[-0.3,0.5],'k','linestyle','--')
    plot([1000,1000],[-0.3,0.5],'k','linestyle','--')
    plot([1250,1250],[-0.3,0.5],'k','linestyle','--')
    hold off;
end

