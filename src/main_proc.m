%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part2Track
% Two-dimensional particle tracking velocimetry (2D-PTV) algorithm for the
% evaluation of double-frame or time-resolved image sequences.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 21.09.17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% include subroutines
addpath(genpath('postprocessing'))
addpath(genpath('tracking'))

% Warnings
warning('off','MATLAB:rankDeficientMatrix')
spmd
warning('off','MATLAB:rankDeficientMatrix')
end

%%
disp('Part2Track.');
%% Measurement Parameters
% folder information
dir_eval = '..\test_cases\double_frame\Lung\';

% define results folder
dir_save = [dir_eval,'results\']; % folder to save results

% load parameter file
run([dir_eval,'parameter.m']);

% plotting options
options.plot_int_results = 1; % Show intermediate results: 1-yes 0-no

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch options.acq_mode
    case 'double_frame'        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Double Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%
        for tt=1:options.n_frames % Initialize Trajectories
            %% import images
            temp_A = double(imread([dir_eval,options.im_pre,num2str(tt,'%.4d'),'A.',options.im_file]));
            temp_B = double(imread([dir_eval,options.im_pre,num2str(tt,'%.4d'),'B.',options.im_file]));

            im_A = temp_A;
            im_B = temp_B;
            
            %% apply image mask
            if isfield(options,'im_roi')
                im_A = im_A.*options.im_roi;
                im_B = im_B.*options.im_roi;
            end
            
            % clear temp variables
            clear temp_A clear temp_B
            
            %% particle detection
            disp('###Perform particle detection###')

            part_A_all = f_detect_particles(double(im_A),...
                                        options.p_size,options.p_int);
            part_B_all = f_detect_particles(double(im_B),...
                                        options.p_size,options.p_int);

            disp(['Particles found in image A: ',num2str(size(part_A_all,1))])
            disp(['Particles found in image B: ',num2str(size(part_B_all,1))])
            
            close(gcf)
            
            %% particle matching
            [results, results_px] = f_track_double_frame(part_A_all,part_B_all,options);
            
            %% save data
            if ~exist(dir_save, 'dir')
                mkdir(dir_save)
            end
            save([dir_save,'results_',num2str(tt),'.mat'],'results');
            save([dir_save,'results_px_',num2str(tt),'.mat'],'results_px');
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Double Frame %%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'time_resolved'        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time Resolved %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Particle Identification
        for tt=1:options.n_frames % Find particle matches in first four frames
            %% import images
            temp_A = double(imread([dir_eval,options.im_pre,num2str(tt,'%.4d'),'.',options.im_file]));

            im_A = temp_A;
            
            %% apply image mask
            if isfield(options,'im_roi')
                im_A = im_A.*options.im_roi;
            end

            %% clear temp variables
            clear temp_A clear temp_B
            
            %% peak detection
            disp('###Perform particle detection###')

            part = f_detect_particles(double(im_A),...
                                        options.p_size,options.p_int);

            %% save particle data
            if ~exist(dir_save, 'dir')
                mkdir(dir_save)
            end
            save([dir_save,'part_',num2str(tt),'.mat'],'part');
        end

        %% Tracking        
        [traj_px] = f_track_time_resolved(dir_save,options);

        %% Save results
        % Apply scaling
        traj_mat = cell2mat(traj_px);
        traj_mat(1:4:end,:) = traj_mat(1:4:end,:)*options.m;
        traj_mat(2:4:end,:) = traj_mat(2:4:end,:)*options.m;
        traj_mat(3:4:end,:) = traj_mat(3:4:end,:)*options.m/options.dt/1000;
        traj_mat(4:4:end,:) = traj_mat(4:4:end,:)*options.m/options.dt/1000;
        
        % Convert matrix to cell array
        traj = mat2cell(traj_mat,[4*ones(size(traj_px,1),1)], [size(traj_mat,2)]);
        
        % Save trajectories
        save([dir_save,'traj_px.mat'],'traj_px');
        save([dir_save,'traj.mat'],'traj');
        
        % Vector fields for all time steps
        traj_mat = cell2mat(traj_px);
        for tt=1:options.n_frames
            results_px = [traj_mat(1:4:end,tt)... % px_x
                          traj_mat(2:4:end,tt)... % px_y
                          traj_mat(3:4:end,tt)... % d_px_x
                          traj_mat(4:4:end,tt)... % d_px_y
                          ones(size(traj_mat(1:4:end,tt),1),1)]; % valid
                      
            results = [results_px(:,1)*options.m... % x
                       results_px(:,2)*options.m... % y
                       results_px(:,3)*options.m/options.dt/1000 ... % u
                       results_px(:,4)*options.m/options.dt/1000 ... % v 
                       results_px(:,5)]; % valid
     
            % Save result matrices for each time step
            save([dir_save,'results_',num2str(tt),'.mat'],'results');
            save([dir_save,'results_px_',num2str(tt),'.mat'],'results_px');
        end
        
        %% Final result
        f_post_plot_trajectories(traj,options);
        xlim([0 options.im_res(2)*options.m])
        ylim([0 options.im_res(1)*options.m])
        set(gca,'ydir','reverse')
        cb = colorbar;
        box on
        set(gca,'FontSize',12);        
        axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time Resolved %%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Warnings
warning('on','MATLAB:rankDeficientMatrix')
spmd
warning('on','MATLAB:rankDeficientMatrix')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import result matrices
for ii=1:options.n_frames
    load([dir_save,'results_',num2str(ii),'.mat']);
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
% Vector Plot
figure('name','Vector Field','NumberTitle','off','Color','w');
quiver(x(1:1:end),y(1:1:end),u(1:1:end),v(1:1:end),1.5,'Color','black')
xlim([0 options.im_res(2)*options.m])
ylim([0 options.im_res(1)*options.m])
set(gca,'ydir','reverse')
axis equal
axis off

%% Check for peak locking
plot_var = u; % Variable to plot
plot_var = plot_var/options.m*options.dt*1000; % Scale to pixel shift
f_post_check_peak_locking(plot_var)

%% Interpolate on grid
grid_win = 8;
iw_mult = 2;

meas_grid = struct('X_vol_min', grid_win*options.m, ...
                   'X_vol_max', options.im_res(2)*options.m, ...
                   'dX_vol', grid_win*options.m,...
                   'Y_vol_min', grid_win*options.m,...
                   'Y_vol_max', options.im_res(1)*options.m,...
                   'dY_vol', grid_win*options.m);
               
[ u_mean, v_mean, u_abs_mean, u_all, v_all ,grid_count,w_all] = ...
    f_post_create_velo_field_gauss_grid_based(x(valid_flag),y(valid_flag),...
                                              u(valid_flag),v(valid_flag), meas_grid, iw_mult);

% Calculate TKE
[ tke, u_fluc_2_mean, v_fluc_2_mean, u_v_fluc_mean] = ...
    f_post_calc_tke(u_mean,v_mean, u_all, v_all,w_all, meas_grid);
%% Show results

%Plot Grid Count
plot_var = grid_count';
if exist('im_roi','var')
	plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
end
figure('name','Number of vectors per grid point','NumberTitle','off','Color','w');
imagesc((plot_var))
colorbar
shading interp
caxis([0 max(plot_var(:))])

% Plot Velocity Field
plot_var = u_abs_mean';
if exist('im_roi','var')
	plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = NaN;
end
figure('name','Mean Velocity Field','NumberTitle','off','Color','w');
imagesc((plot_var))
caxis([0 max(plot_var(:))])
colorbar
shading interp
axis off

% Plot Velocity Field
plot_var = tke';

if exist('im_roi','var')
 plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
end
figure('name','TKE Field','NumberTitle','off','Color','w');
imagesc((plot_var))
colorbar
shading interp
caxis([0 max(plot_var(:))])
colorbar
shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%