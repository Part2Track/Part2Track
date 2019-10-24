function [ u_grid_mean, v_grid_mean, u_abs_grid_mean, u_grid_all, v_grid_all, grid_count] = f_post_create_velo_field_polynomial_grid_based(x,y,u,v,meas_grid,r_o_a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of scattered PTV data on regular grid by looping through 
% all grid points (collecting-vectors-style).
% Use polynomial interpolation to find velocity value at center spot based
% on:
% N. Agüera, G. Cafiero, T. Astarita, S. Discetti: "Ensemble 3D PTV for 
% high resolution turbulent statistics"
% Measurement Science and Technology, 2016, 27, 124011
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
% origin: Thomas Janke / 26.07.19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nx = meas_grid.X_vol_min:meas_grid.dX_vol:meas_grid.X_vol_max;
ny = meas_grid.Y_vol_min:meas_grid.dY_vol:meas_grid.Y_vol_max;
u_grid_mean =zeros(size(nx,2),size(ny,2));
v_grid_mean =zeros(size(nx,2),size(ny,2));
u_abs_grid_mean =zeros(size(nx,2),size(ny,2));
grid_count =zeros(size(nx,2),size(ny,2));
radius = sqrt((meas_grid.dX_vol)^2+(meas_grid.dY_vol)^2)*r_o_a;
u_grid_all = cell(size(nx,2),size(ny,2));
v_grid_all = cell(size(nx,2),size(ny,2));
n_nx = size(nx,2);
n_ny = size(ny,2);

parfor ind_x =1:n_nx
    find_x = (x>(nx(ind_x)-radius)) & (x<(nx(ind_x)+radius));
    x_temp_x = x(find_x);
    y_temp_x = y(find_x);
    u_temp_x = u(find_x);
    v_temp_x = v(find_x);
    
    for ind_y =1:n_ny
        find_y = (y_temp_x>(ny(ind_y)-radius)) & (y_temp_x<(ny(ind_y)+radius));
        x_temp = x_temp_x(find_y);
        y_temp = y_temp_x(find_y);
        u_temp = u_temp_x(find_y);
        v_temp = v_temp_x(find_y);

        x_fit = x_temp-nx(ind_x);
        y_fit = y_temp-ny(ind_y);
        
        if length(x_fit) > 6
            sf_u = fit([x_fit, y_fit],u_temp,'poly22');
            u_temp_fit = feval(sf_u,[0 0]);
            sf_v = fit([x_fit, y_fit],v_temp,'poly22');
            v_temp_fit = feval(sf_v,[0 0]);

            u_abs_temp = sqrt(u_temp_fit.^2+v_temp_fit.^2);

            u_grid_mean(ind_x,ind_y) = u_temp_fit;
            v_grid_mean(ind_x,ind_y) = v_temp_fit;
            u_abs_grid_mean(ind_x,ind_y) = u_abs_temp;

            u_grid_all{ind_x,ind_y} = u_temp;
            v_grid_all{ind_x,ind_y} = v_temp;

            grid_count(ind_x,ind_y) = size(u_temp,1);
        end
        
    end
end

u_grid_mean(u_grid_mean==0) = NaN;
v_grid_mean(v_grid_mean==0) = NaN;
u_abs_grid_mean(u_abs_grid_mean==0) = NaN;
grid_count(grid_count==0) = NaN;

end

