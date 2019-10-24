function [ valid_flag ] = f_post_outlier_median( x,y,u,v,n_neighbours,thr,noise )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal outlier detection method for scattered data based on:
% J. Westerweel, F. Scarano: "Univsersal outlier detection for PIV data"
% Experiments in Fluids, 2005, 39:1096 
%
% Parallelized version of this algorithm adopted from:
% M. Patel, S.E. Leggett, A.K. Landauer, I.Y. Wong, C. Franck:
% "Rapid, topology-based particle tracking for high-resolution measurements
% of large complex 3D motion fields"
% Scientific Reports, 2018, 8, Article number:5581
% Github: https://github.com/FranckLab/T-PT
% under MIT License (see below)
% Copyright (c) 2018 Franck Lab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      x             - x-coordinate of vectors [n x 1]
%   ------
%               y             - y-coordinate of vectors [n x 1]
% 
%               u             - x-velocity-component of vector [n x 1]
% 
%               v             - y-velocity-component of vector [n x 1]
%
%               n_neighbours  - number of neighbouring particles [integer]
%
%               thr           - fluctuation threshold for 
%                               outlier detection [integer]
%                               |
%                               ---> recommended: thr = 2
%
%               noise          - noise level of measurement [integer]
%                                |
%                                ---> recoomended: noise = 0.1
%
%
%   Output:     valid_flag     - logical vector of outlier detection choice
%   -------                         |
%                                   |--> 1: valid vector
%                                   |
%                                   |--> 0: invalid vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 17.11.17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Serial 
% % initialize variables
% norm_fluc = zeros(size(x,1),2); 
% 
% % loop over velocity vectors
% for ii=1:size(x,1)        
%     
%     % find neighbouring vectors
%     temp_dist = sqrt((x(:)-x(ii)).^2+(y(:)-y(ii)).^2);
%     [~, index_neighbours] = sort(temp_dist);
%     % discard vectors too far away
%     index_neighbours(n_neighbours+2:end) = 0;
%     index_neighbours(index_neighbours == 0) = [];
%     % discard current vector from list
%     index_neighbours(1) = [];
%     
%     % loop over velocity components
%     for cc = 1:2
%         if cc == 1
%             vel = u;
%         else
%             vel = v;
%         end
%      med_vel = median(vel(index_neighbours)); %median of neighbouring vectors
%      vel_fluc = vel(ii)-med_vel; %fluctuation of current vector with respect to median
%      vel_res = vel(index_neighbours)-med_vel; %residual: fluctuations of neighbouring vectors with respect to median
%      med_res = median(abs(vel_res)); %median of absolute residuals
%      norm_fluc(ii,cc) = abs(vel_fluc/(med_res+noise)); %normalised fluctuations
%     end
% end    
% 
% valid_flag = sqrt(norm_fluc(:,1).^2 + norm_fluc(:,2).^2) < thr; %detection criterion


%% Vectorized
% 
idx = knnsearch([x y],[x y],'k',n_neighbours); % Find n_neighbours particles for every displacement
idx = idx(:,2:end); % Remove first point. It's the self point

un1 = zeros(size(idx)); % Initialize neighbour matrix
vn2 = un1; % Initialize neighbour matrix

% median of neighbouring vectors
un1(:) = u(idx(:)); med_vel_x = median(un1,2);
vn2(:) = v(idx(:)); med_vel_y = median(vn2,2);

% residual: fluctuations of neighbouring vectors with respect to median
vel_res_x = abs(bsxfun(@plus,un1,-med_vel_x)); 
vel_res_y = abs(bsxfun(@plus,vn2,-med_vel_y)); 

% median of absolute residuals
med_res_x = median(vel_res_x,2);
med_res_y = median(vel_res_y,2);

% fluctuation of investigated vector with respect to median
vel_fluc_x = abs(bsxfun(@plus,u,-med_vel_x)); 
vel_fluc_y = abs(bsxfun(@plus,v,-med_vel_y)); 

% normalised fluctuations
norm_fluc_x = vel_fluc_x./(med_res_x+noise);
norm_fluc_y = vel_fluc_y./(med_res_y+noise);

% detection criterion
valid_flag = sqrt(norm_fluc_x.^2+norm_fluc_y.^2) < thr;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
% Copyright (c) 2018 Franck Lab
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
