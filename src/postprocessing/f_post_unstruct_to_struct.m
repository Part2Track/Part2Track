function [] = f_post_unstruct_to_struct(dir_save,parameter,sample_mode,bin_method,grid_win,meas_grid,iw_mult)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation of scattered PTV data on regular grid for all timesteps.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      dir_save         - path to scattered results
%   ------
%               parameter        - parameter structure
% 
%               acq_mode         - image acquisition mode
%
%               grid_win         - grid spacing
%
%               meas_grid        - structure of grid information
%                   |
%                   |----> 'X_vol_min' 
%                   |----> 'X_vol_max'
%                   |----> 'dX_vol'
%                   |----> 'Y_vol_min'
%                   |----> 'Y_vol_max'
%                   |----> 'dY_vol'
%
%               iw_mult          - multiplier for IW size
%                                 (IW=grid_win*iw_mult)
%
%   Output:     Structured velocity fields are saved in the subdirectory 
%   -------     'structured'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 29.05.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

n_frames = parameter.n_frames; % number of double-frames 
acq_mode = parameter.acq_mode;
time_avg = 0;

if isfield(parameter,'im_roi')
    im_roi = parameter.im_roi; % import image mask, if defined
end

switch sample_mode
    case 'single'
        parfor ii=1:n_frames
            if strcmp(acq_mode,'double_frame') == 1
                range_load = [ii-time_avg:1:ii+time_avg];
                if range_load(1) < 1
                    range_load = range_load + abs(range_load(1)) + 1;
                elseif range_load(end) > n_frames
                    range_load = range_load - (range_load(end)-n_frames);
                end
                
                for jj=1:size(range_load,2)       
                    range_load(jj)
                    results=load([dir_save,'results_',num2str(range_load(jj)),'.mat']);
                    results = results.results;
                    if jj==1
                        x_scat = results(:,1);
                        y_scat = results(:,2);
                        u_scat = results(:,3);
                        v_scat = results(:,4);
                        valid_flag = logical(results(:,5));
                    else
                        x_scat = [x_scat; results(:,1)];
                        y_scat = [y_scat; results(:,2)];
                        u_scat = [u_scat; results(:,3)];
                        v_scat = [v_scat; results(:,4)];
                        valid_flag = [valid_flag; logical(results(:,5))];
                    end
                end
            else
                results=load([dir_save,'results_',num2str(ii),'.mat']);
                results = results.results;

                x_scat = results(:,1);
                y_scat = results(:,2);
                u_scat = results(:,3);
                v_scat = results(:,4);
                valid_flag = logical(results(:,5));                
            end

           [ u_temp, v_temp, u_abs_temp, ~, ~ ,~] = f_post_create_velo_field_gauss_grid_based( x_scat(valid_flag),y_scat(valid_flag),u_scat(valid_flag),v_scat(valid_flag), meas_grid, iw_mult);

            u_temp = u_temp'.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2));
            v_temp = v_temp'.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2));

            u_abs_temp = sqrt(u_temp.^2+v_temp.^2);

            u(:,:,ii) = u_temp;
            v(:,:,ii) = v_temp;
            u_abs(:,:,ii) = u_abs_temp;
            grid_temp(:,:,ii) = ~isnan(u_temp);
        end
        
        u_mean = nanmean(u,3);
        v_mean = nanmean(v,3);
        u_abs_mean = nanmean(u_abs,3);

        u_std = nanstd(u,1,3);
        v_std = nanstd(v,1,3);
        u_abs_std = nanstd(u_abs,1,3);

        [tke, u_fluc_2_mean, v_fluc_2_mean, u_v_fluc_mean] = f_post_calc_tke(u_mean,v_mean, u, v, [], meas_grid);
        
        grid_count = nansum(double(grid_temp),3);


        %% save data
        if ~exist([dir_save,'structured\'], 'dir')
            mkdir([dir_save,'structured\'])
        end

        save([dir_save,'structured\grid_count.mat'],'grid_count');

        save([dir_save,'structured\u.mat'],'u');
        save([dir_save,'structured\v.mat'],'v');
        save([dir_save,'structured\u_abs.mat'],'u_abs');

        save([dir_save,'structured\u_mean.mat'],'u_mean');
        save([dir_save,'structured\v_mean.mat'],'v_mean');
        save([dir_save,'structured\u_abs_mean.mat'],'u_abs_mean');

        save([dir_save,'structured\u_std.mat'],'u_std');
        save([dir_save,'structured\v_std.mat'],'v_std');
        save([dir_save,'structured\u_abs_std.mat'],'u_abs_std');

        save([dir_save,'structured\tke.mat'],'tke');
        save([dir_save,'structured\u_fluc_2_mean.mat'],'u_fluc_2_mean');
        save([dir_save,'structured\v_fluc_2_mean.mat'],'v_fluc_2_mean');
        save([dir_save,'structured\u_v_fluc_mean.mat'],'u_v_fluc_mean');
        
    case 'ensemble'
        for ii=1:n_frames
            results=load([dir_save,'results_',num2str(ii),'.mat']);
            results = results.results;
            if ii == 1
                x = results(:,1);
                y = results(:,2);
                u = results(:,3);
                v = results(:,4);
                valid_flag = logical(results(:,5));
            else
                x = [x; results(:,1)];
                y = [y; results(:,2)];
                u = [u; results(:,3)];
                v = [v; results(:,4)];
                valid_flag = [valid_flag; logical(results(:,5))];
            end
        end
        
        % Calculate mean velocity fields
        switch bin_method
            case 'Gaussian'
                [u_mean, v_mean, u_abs_mean,u, v ,grid_count, w_all] = ...
                    f_post_create_velo_field_gauss_grid_based( x(valid_flag),...
                                                               y(valid_flag),...
                                                               u(valid_flag),...
                                                               v(valid_flag),...
                                                               meas_grid, iw_mult);
            case 'polynomial'               
                [u_mean, v_mean, u_abs_mean,u, v ,grid_count] = ...
                   f_post_create_velo_field_polynomial_grid_based( x(valid_flag),...
                                                                   y(valid_flag),...
                                                                   u(valid_flag),...
                                                                   v(valid_flag), ...
                                                                   meas_grid, iw_mult);
                                                              
        end

        % Apply image mask
        u_mean = u_mean.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';        
        v_mean = v_mean.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';      
        u_abs_mean = u_abs_mean.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';
              
        grid_count = grid_count.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';
        
        % Calculate velocity fluctuation fields
        switch bin_method
            case 'Gaussian'
                [tke, u_fluc_2_mean, v_fluc_2_mean, u_v_fluc_mean] = ...
                    f_post_calc_tke(u_mean,v_mean, u, v, w_all, meas_grid);
            case 'polynomial'
                [tke, u_fluc_2_mean, v_fluc_2_mean, u_v_fluc_mean] = ...
                    f_post_calc_tke(u_mean,v_mean, u, v, [], meas_grid);
        end
                 
        % Apply image mask
        tke = tke.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';      
        u_fluc_2_mean = u_fluc_2_mean.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';      
        v_fluc_2_mean = v_fluc_2_mean.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))';      
        u_v_fluc_mean = u_v_fluc_mean.*im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2))'; 
        
        % Transpose matrices
        grid_count = grid_count';

        u_mean = u_mean';
        v_mean = v_mean';
        u_abs_mean = u_abs_mean';

        tke = tke';
        u_fluc_2_mean = u_fluc_2_mean';
        v_fluc_2_mean = v_fluc_2_mean';
        u_v_fluc_mean = u_v_fluc_mean';

        % Save data
        if ~exist([dir_save,'structured\'], 'dir')
            mkdir([dir_save,'structured\'])
        end
        save([dir_save,'structured\grid_count.mat'],'grid_count');

        save([dir_save,'structured\u_mean.mat'],'u_mean');
        save([dir_save,'structured\v_mean.mat'],'v_mean');
        save([dir_save,'structured\u_abs_mean.mat'],'u_abs_mean');

        save([dir_save,'structured\tke.mat'],'tke');
        save([dir_save,'structured\u_fluc_2_mean.mat'],'u_fluc_2_mean');
        save([dir_save,'structured\v_fluc_2_mean.mat'],'v_fluc_2_mean');
        save([dir_save,'structured\u_v_fluc_mean.mat'],'u_v_fluc_mean');
end

end

