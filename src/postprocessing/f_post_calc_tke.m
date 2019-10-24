function [ tke, u_fluc_2_mean, v_fluc_2_mean, u_v_fluc_mean ] = f_post_calc_tke(u_grid_mean,v_grid_mean, u_grid_all, v_grid_all, w_all,meas_grid )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Turbulent Kinetic Energy and 
% Mean Velocity Fluctuations <u'u'>, <v'v'> & <u'v'>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:    | u_grid_mean        - mean u-velocity field 
%   ------    |
%             | v_grid_mean        - mean v-velocity field
%             |
%             | u_grid_all         - instantaneous u-velocity fields
%             |
%             | v_grid_all         - instantaneous u-velocity fields
%             |
%             | w_all              - weightings for averaging
%             |
%             |-----> all matrices from:
%                     f_post_create_velo_field_gauss_grid_based
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
%   Output:  tke       - Turbulent Kinetic Energy [size(nx,2) x size(ny,2)]
%   -------
%            u_fluc_2_mean   - <u'u'> [size(nx,2) x size(ny,2)]
%
%            v_fluc_2_mean   - <v'v'> [size(nx,2) x size(ny,2)]
%
%            u_v_fluc_mean   - <u'v'> [size(nx,2) x size(ny,2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 28.03.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nx = meas_grid.X_vol_min:meas_grid.dX_vol:meas_grid.X_vol_max;
ny = meas_grid.Y_vol_min:meas_grid.dY_vol:meas_grid.Y_vol_max;

% Calculate Mean Fluctuations and Turbulent Kinetic Energy
if isempty(w_all) == 0
    
u_fluc_2 = cell(size(nx,2),size(ny,2));
v_fluc_2 = cell(size(nx,2),size(ny,2));
u_v_fluc = cell(size(nx,2),size(ny,2));

    for ind_x =1:size(nx,2)   
        for ind_y =1:size(ny,2)
            if ~isnan(u_grid_mean(ind_x,ind_y))
                u_fluc_2{ind_x,ind_y} = (repmat(u_grid_mean(ind_x,ind_y),size(u_grid_all{ind_x, ind_y},1),1)-...
                                                u_grid_all{ind_x, ind_y}).^2.*w_all{ind_x, ind_y};
                v_fluc_2{ind_x,ind_y} = (repmat(v_grid_mean(ind_x,ind_y),size(v_grid_all{ind_x, ind_y},1),1)-...
                                                v_grid_all{ind_x, ind_y}).^2.*w_all{ind_x, ind_y};
                u_v_fluc{ind_x,ind_y} = (repmat(u_grid_mean(ind_x,ind_y),size(u_grid_all{ind_x, ind_y},1),1)-...
                                                u_grid_all{ind_x, ind_y}).*...
                                             (repmat(v_grid_mean(ind_x,ind_y),size(v_grid_all{ind_x, ind_y},1),1)-...
                                                v_grid_all{ind_x, ind_y}.*w_all{ind_x, ind_y});
            else
                u_fluc_2{ind_x, ind_y} = NaN;
                v_fluc_2{ind_x, ind_y} = NaN;
                u_v_fluc{ind_x,ind_y} = NaN;
            end
        end
    end

        sum_w = cellfun(@nansum, w_all);

        u_fluc_2_w = cellfun(@nansum, u_fluc_2);
        v_fluc_2_w = cellfun(@nansum, v_fluc_2);
        u_v_fluc_w = cellfun(@nansum, u_v_fluc);

        u_fluc_2_mean = u_fluc_2_w./sum_w;
        v_fluc_2_mean = v_fluc_2_w./sum_w;
        u_v_fluc_mean = u_v_fluc_w./sum_w;
        tke = 3/4*(u_fluc_2_mean+v_fluc_2_mean);

else
    
    if iscell(u_grid_all) == 0
        u_fluc_2_mean = nanmean((repmat(u_grid_mean,1,1,size(u_grid_all,3))-u_grid_all).^2,3);
        v_fluc_2_mean = nanmean((repmat(v_grid_mean,1,1,size(v_grid_all,3))-v_grid_all).^2,3);
        u_v_fluc_mean = nanmean((repmat(v_grid_mean,1,1,size(v_grid_all,3))-v_grid_all).*(repmat(u_grid_mean,1,1,size(u_grid_all,3))-v_grid_all),3);
        tke = 3/4*(u_fluc_2_mean+v_fluc_2_mean);
    else

        u_fluc_2 = cell(size(nx,2),size(ny,2));
        v_fluc_2 = cell(size(nx,2),size(ny,2));
        u_v_fluc = cell(size(nx,2),size(ny,2));

        for ind_x =1:size(nx,2)   
            for ind_y =1:size(ny,2)
                if ~isnan(u_grid_mean(ind_x,ind_y))
                    u_fluc_2{ind_x,ind_y} = (repmat(u_grid_mean(ind_x,ind_y),size(u_grid_all{ind_x, ind_y},1),1)-...
                                                    u_grid_all{ind_x, ind_y}).^2;
                    v_fluc_2{ind_x,ind_y} = (repmat(v_grid_mean(ind_x,ind_y),size(v_grid_all{ind_x, ind_y},1),1)-...
                                                    v_grid_all{ind_x, ind_y}).^2;
                    u_v_fluc{ind_x,ind_y} = (repmat(u_grid_mean(ind_x,ind_y),size(u_grid_all{ind_x, ind_y},1),1)-...
                                                    u_grid_all{ind_x, ind_y}).*...
                                                 (repmat(v_grid_mean(ind_x,ind_y),size(v_grid_all{ind_x, ind_y},1),1)-...
                                                    v_grid_all{ind_x, ind_y});
                else
                    u_fluc_2{ind_x, ind_y} = NaN;
                    v_fluc_2{ind_x, ind_y} = NaN;
                    u_v_fluc{ind_x,ind_y} = NaN;
                end
            end
        end

            u_fluc_2_mean = cellfun(@nanmean, u_fluc_2);
            v_fluc_2_mean = cellfun(@nanmean, v_fluc_2);
            u_v_fluc_mean = cellfun(@nanmean, u_v_fluc);

            tke = 3/4*(u_fluc_2_mean+v_fluc_2_mean);
   end
end
end

