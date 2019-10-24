function [ u_grid_mean, v_grid_mean, u_abs_grid_mean, u_grid_all, v_grid_all, grid_count, w_all] = f_post_create_velo_field_gauss_grid_based(x,y,u,v,meas_grid,r_o_a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of scattered PTV data on regular grid by looping through 
% all grid points (collecting-vectors-style).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      x         - x-position of vector [n x 1]
%   ------
%               y         - y-position of vector [n x 1]
% 
%               u         - u-velocity of vector [n x 1]
%
%               v         - v-velocity of vector [n x 1]
%
%               meas_grid - structure of grid information
%                   |
%                   |----> 'X_vol_min' 
%                   |----> 'X_vol_max'
%                   |----> 'dX_vol'
%                   |----> 'Y_vol_min'
%                   |----> 'Y_vol_max'
%                   |----> 'dY_vol'
%
%               r_o_a      - multiplier for radius of averaging area
%
%   Output:     u_grid     - u-velocity on grid [size(nx,2) x size(ny,2)]
%   -------
%               v_grid     - v-velocity on grid [size(nx,2) x size(ny,2)]
%
%               grid_count - counts of samples within grid-point
%                            [size(nx,2) x size(ny,2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 21.09.17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nx = meas_grid.X_vol_min:meas_grid.dX_vol:meas_grid.X_vol_max;
ny = meas_grid.Y_vol_min:meas_grid.dY_vol:meas_grid.Y_vol_max;
sum_w =zeros(size(nx,2),size(ny,2));
u_grid_w =zeros(size(nx,2),size(ny,2));
v_grid_w =zeros(size(nx,2),size(ny,2));
u_abs_grid_w =zeros(size(nx,2),size(ny,2));
grid_count =zeros(size(nx,2),size(ny,2));
radius = sqrt((meas_grid.dX_vol)^2+(meas_grid.dY_vol)^2)*r_o_a;
sigma = 0.5*radius/1;
u_grid_all = cell(size(nx,2),size(ny,2));
v_grid_all = cell(size(nx,2),size(ny,2));
w_all = cell(size(nx,2),size(ny,2));

for ind_x =1:size(nx,2)
    find_x = (x>(nx(ind_x)-radius)) & (x<(nx(ind_x)+radius));
    x_temp_x = x(find_x);
    y_temp_x = y(find_x);
    u_temp_x = u(find_x);
    v_temp_x = v(find_x);
    
    for ind_y =1:size(ny,2)
        find_y = (y_temp_x>(ny(ind_y)-radius)) & (y_temp_x<(ny(ind_y)+radius));
        x_temp = x_temp_x(find_y);
        y_temp = y_temp_x(find_y);
        u_temp = u_temp_x(find_y);
        v_temp = v_temp_x(find_y);
        
        u_abs_temp = sqrt(u_temp.^2+v_temp.^2);
        dist_vect = sqrt((x_temp-nx(ind_x)).^2 + (y_temp-ny(ind_y)).^2);
        
        w = 1/(sigma*sqrt(2*pi))*exp(-dist_vect.^2/(2*(sigma)^2));

        sum_w(ind_x,ind_y) = nansum(w);

        u_grid_w(ind_x,ind_y) = nansum(u_temp.*w);
        v_grid_w(ind_x,ind_y) = nansum(v_temp.*w);
        u_abs_grid_w(ind_x,ind_y) = nansum(u_abs_temp.*w);
        
        u_grid_all{ind_x,ind_y} = u_temp;
        v_grid_all{ind_x,ind_y} = v_temp;
        
        grid_count(ind_x,ind_y) = size(w,1);
        
        w_all{ind_x,ind_y} = w;
    end
end

% Normalize with sum of weights
u_grid_mean = u_grid_w./sum_w;
u_grid_mean(u_grid_mean==0) = NaN;
v_grid_mean = v_grid_w./sum_w;
v_grid_mean(v_grid_mean==0) = NaN;
u_abs_grid_mean = u_abs_grid_w./sum_w;
u_abs_grid_mean(u_abs_grid_mean==0) = NaN;
grid_count(grid_count==0) = NaN;

end

