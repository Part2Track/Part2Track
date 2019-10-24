function [ u_grid, v_grid, grid_count] = f_post_create_velo_field_gauss( x,y,u,v, meas_grid )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of scattered PTV data on regular grid by looping through 
% all vectors (pick-and-place-style).
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
grid_count =zeros(size(nx,2),size(ny,2));

radius = sqrt((meas_grid.dX_vol)^2+(meas_grid.dY_vol)^2);

for ii = 1:size(u,1)            
    ind_x =find(( abs(nx(:)-x(ii))== min( abs(nx(:)-x(ii)))));
    if size(ind_x,1) > 1
        ind_x(2) =[];
    end
    
    ind_y =find(( abs(ny(:)-y(ii))== min( abs(ny(:)-y(ii)))));
    if size(ind_y,1) > 1
        ind_y(2) =[];
    end

    center_point = [meas_grid.X_vol_min+(ind_x-1)*meas_grid.dX_vol meas_grid.Y_vol_min+(ind_y-1)*meas_grid.dY_vol];
    dist_vect = sqrt( (x(ii)-center_point(1))^2 + (y(ii)-center_point(2))^2 );
    if dist_vect < radius
%     w = exp(-dist_vect^2/sqrt(-radius^2/(2*log(0.01))));
    w = exp(-(dist_vect/radius*2)^2/sqrt(-(radius/radius)^2/(2*log(0.01))));
    sum_w(ind_x,ind_y) = sum_w(ind_x,ind_y)+w;

    u_grid_w(ind_x,ind_y) = u_grid_w(ind_x,ind_y)+u(ii)*w;
    v_grid_w(ind_x,ind_y) = v_grid_w(ind_x,ind_y)+v(ii)*w;


    grid_count(ind_x,ind_y) = grid_count(ind_x,ind_y) +1;
    end     
end
u_grid = u_grid_w./sum_w;
u_grid(u_grid==0) = NaN;
v_grid = v_grid_w./sum_w;
v_grid(v_grid==0) = NaN;
grid_count(grid_count==0) = NaN;

end

