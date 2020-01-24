%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main post-processing routine for Part2Track.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 29.05.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 

% include subroutines
addpath(genpath('src\postprocessing'))
addpath(genpath('src\tracking'))

%%
disp('Part2Track Post Processing.');

%% Parameters
% folder information
dir_eval = '..\test_cases\double_frame\Lung\';
dir_save = [dir_eval,'results\']; % folder to save results

% image parameters
run([dir_eval,'parameter.m']);

im_res = options.im_res; % camera resolution
n_frames = options.n_frames; % number of double-frames 
dt = options.dt; % time separation of image pairs [s]
m = options.m; % mapping scale [mm/px]
im_roi = options.im_roi;

% post processing parameter
create_struct = true ; % create structured velocity field (true/false)
sample_mode = 'ensemble'; % 'single' or 'ensemble'
bin_method = 'Gaussian'; % 'Gaussian' or 'polynomial' (only for 'ensemble')
grid_win = 8; % grid-spacing
iw_mult = 2; % interrogation window multiplier (IW = grid_win*iw_mult)
meas_grid = struct('X_vol_min', grid_win*m, 'X_vol_max', im_res(2)*m, 'dX_vol', grid_win*m, 'Y_vol_min', grid_win*m, 'Y_vol_max', im_res(1)*m, 'dY_vol', grid_win*m); % create structured grid

%% Convert unstructured velocity fields to structured velocity fields
if create_struct == true
    f_post_unstruct_to_struct(dir_save,options,sample_mode,bin_method,grid_win,meas_grid,iw_mult);
end
%% Load structured data
switch sample_mode
    case 'single'
        load([dir_save,'structured\grid_count']); 
        
        load([dir_save,'structured\u']); 
        load([dir_save,'structured\v']); 
        load([dir_save,'structured\u_abs']); 

        load([dir_save,'structured\u_mean']); 
        load([dir_save,'structured\v_mean']); 
        load([dir_save,'structured\u_abs_mean']); 

        load([dir_save,'structured\u_std']); 
        load([dir_save,'structured\v_std']); 
        load([dir_save,'structured\u_abs_std']); 

        load([dir_save,'structured\tke']); 
        load([dir_save,'structured\u_fluc_2_mean']); 
        load([dir_save,'structured\v_fluc_2_mean']); 
        load([dir_save,'structured\u_v_fluc_mean']); 
        
    case 'ensemble'
        load([dir_save,'structured\grid_count']); 
        
        load([dir_save,'structured\u_mean']); 
        load([dir_save,'structured\v_mean']); 
        load([dir_save,'structured\u_abs_mean']); 

        load([dir_save,'structured\tke']); 
        load([dir_save,'structured\u_fluc_2_mean']); 
        load([dir_save,'structured\v_fluc_2_mean']); 
        load([dir_save,'structured\u_v_fluc_mean']); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collection of various post-processing routines
% Choose what you want to process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Smooth Kernel
smooth_flag = false;
filtWidth = 2;
filtSigma = 5;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
%% Grid Count
figure('name','Number of vectors per grid point','NumberTitle','off');
plot_var = grid_count;
if exist('im_roi','var')
	plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
end
imagesc((plot_var))
colorbar
shading interp
caxis([0 max(plot_var(:))])
set(gcf,'color','w');
axis off
%% Magnitude
plot_var = u_abs_mean;
if exist('im_roi','var')
 plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = NaN;
end
if smooth_flag == true
 plot_var = nanconv(plot_var,imageFilter, 'nanout');
end
figure('name','Mean Velocity Field','NumberTitle','off');   
set(gcf,'color','w');
imagesc((plot_var))
colorbar
shading interp
caxis([0 max(plot_var(:))])
axis off
set(gca,'YDir','reverse')

%% Turbulent Kinetic Energy
plot_var = tke;
if smooth_flag == true
 plot_var = nanconv(plot_var,imageFilter, 'nanout');
end
if exist('im_roi','var')
 plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
end

figure('name','Turbulent Kinetic Energy','NumberTitle','off');
imagesc((plot_var))
colorbar
shading interp
caxis([0 max(plot_var(:))])
set(gcf,'color','w');
axis off

%% Check for peak locking
plot_var = u_abs_mean; % Variable to plot
plot_var = plot_var/m*dt*1000; % Scale to pixel shift
f_post_check_peak_locking(plot_var)

%% Scatter Plot of Velocity
plot_u = u_mean(:,:);
plot_v = v_mean(:,:);
f_post_scatter_u_v(plot_u,plot_v);

%% Calculate Vortex Critera
[Delta,Lambda_ci,Q,Lambda_2] = f_post_vortex_criteria(u,v,grid_win*m,'true');

%% Simple Animation for Vector Fields (only usable, when sample_mode = 'single')
% Variable to plot
plot_var = Q;

% Create new video object & set options
mov_name = ['simple_animation.avi'];
vid = VideoWriter(mov_name);
vid.FrameRate = 30;
open(vid)

% Create figure & set options
f=figure('name','Mean Velocity Field Animation','NumberTitle','off');  
set(f, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
set(gcf,'color','black'); 

for ii=1:size(plot_var,3) % Loop over all time steps
    % Load data from matrix
    plot_mat = plot_var(:,:,ii);

    % Apply mask, if it exists
    if exist('im_roi','var')
     plot_mat(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
    end

    % Fill gaps and set background to invisible
    plot_mat(plot_mat == 0 ) = -NaN;
    plot_mat(plot_mat < 0 ) = 0;
    
    Vq = plot_mat;

    % Additional smoothing with Gauss filter
    if smooth_flag == true
     Vq = nanconv(Vq,imageFilter, 'nanout');
    end

    if ii==1 % Initialize Figure
        % Plot Field
        h=pcolor((Vq));        
        shading interp

        % Rotating for correct display
        set(gca,'xdir','reverse')
        rotate(h,[0 0 1],-180)

        % Create colorbar
        cb=colorbar;
        set(cb,'TickLabelInterpreter','latex','Color','white','Location','EastOutside','FontSize',24)
        t_cb = get(cb,'title');
        set(t_cb,'String','$|\underline{u}|$ [m/s]','interp','latex','Color','white','FontSize',24)
        
        % Set colorbar limits
        caxis([0 0.05])
        
        % Adjust axis appearence
        axis off
        axis equal
    else % Update matrix with new values from current time step
        set(h, 'CData', Vq);
    end

    % Grab current frame and append to video
    frame = getframe(gcf);
    writeVideo(vid,frame)
end
close(vid) % Close video object to finish saving

%% Simple Animations for Trajectories (only usable for 'time_resolved')
% traj_mat = traj_test;
plot_length = 7; % Length of trajectories
mov_name = ['simple_traj_animation.avi'];

% Create Video Object
vid = VideoWriter(mov_name);
vid.FrameRate = 5;
vid.Quality = 90;
open(vid)

% Create Figure
f=figure('name','Trajectories','NumberTitle','off','Color','white');
set(gcf, 'Units', 'Normalized','Position', [0.1 0.1 0.4 0.6]);

% Create first trajectories
traj_mat = cell2mat(traj);
X_tr_plot = traj_mat(1:4:end,1:plot_length)';
X_tr_plot(:,sum(isnan(X_tr_plot)) == 7) = [];
Y_tr_plot = traj_mat(2:4:end,1:plot_length)';
Y_tr_plot(:,sum(isnan(Y_tr_plot)) == 7) = [];
disp_plot = sqrt(traj_mat(3:4:end,1:plot_length).^2+traj_mat(4:4:end,1:plot_length).^2)';
disp_plot(:,sum(isnan(disp_plot)) == 7) = [];

h=patch([X_tr_plot ; nan(1,size(X_tr_plot,2))],[Y_tr_plot; nan(1,size(X_tr_plot,2))],[(disp_plot); nan(1,size(X_tr_plot,2))],'EdgeColor','interp','FaceColor','none','LineWidth',1.75);
xlim([0 parameter.im_res(2)*parameter.m])
ylim([0 parameter.im_res(1)*parameter.m])
axis equal

% Set colors
caxis([0 4.5])

% Create colorbar
cb=colorbar;
set(cb,'TickLabelInterpreter','latex','Color','black','Location','West','FontSize',18)
ylabel(cb,'$\delta x$ [px]','interp','latex','Color','black','FontSize',18)
% t_cb = get(cb,'title');
% set(t_cb,'String','$\delta x$ [px]','interp','latex','Color','black','FontSize',18)
set(cb, 'YAxisLocation','left')

% Change appearance
box off
axis off
set(gca,'FontSize',12); 
set(gca,'ydir','reverse')
set(gca, 'LooseInset', get(gca,'TightInset'))

% Get Frame for Video Object
frame = getframe(gcf);
writeVideo(vid,frame)

for ii=1:parameter.n_frames-plot_length-1 % Loop over time and update plot
        
    % Load new values
    X_tr_plot = traj_mat(1:4:end,1+ii:1+ii+plot_length)';
    X_tr_plot(:,sum(isnan(X_tr_plot)) == 7) = [];
    X_tr_plot = [X_tr_plot ; nan(1,size(X_tr_plot,2))];
    Y_tr_plot = traj_mat(2:4:end,1+ii:1+ii+plot_length)';
    Y_tr_plot(:,sum(isnan(Y_tr_plot)) == 7) = [];
    Y_tr_plot = [Y_tr_plot ; nan(1,size(X_tr_plot,2))];
    disp_plot = sqrt(traj_mat(3:4:end,1+ii:1+ii+plot_length).^2+traj_mat(4:4:end,1+ii:1+ii+plot_length).^2)';
    disp_plot(:,sum(isnan(disp_plot)) == 7) = [];
    disp_plot = [(disp_plot); nan(1,size(X_tr_plot,2))];
    
    % Update plot
    set(h,'XData',X_tr_plot,'YData',Y_tr_plot,'CData',disp_plot)
    drawnow

    % Get Frame for Video Object
    frame = getframe(gcf);
    writeVideo(vid,frame)
end
    
close(vid) % Close video object to finish saving

%% POD decomposition (only usable for 'time_resolved' and sample_mode = 'single')
var_pod = u_abs;
[Psy, Sigma, Phi] = f_post_POD_decompose(var_pod);

%% Show energy content of all modes
figure('name','Turbulent Kinetic Energy','NumberTitle','off');
set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.5]);
subplot(1,2,1)
plot(Sigma./sum(Sigma)*100)
subplot(1,2,2)
plot(cumsum(Sigma./sum(Sigma)*100))
%% Show mode
i_mode = 2; % Select mode
plot_var = Sigma(i_mode)*Psy(:,i_mode);
plot_var = reshape(plot_var,size(var_pod,1),size(var_pod,2));
if smooth_flag == true
 plot_var = nanconv(plot_var,imageFilter, 'nanout');
end
if exist('im_roi','var')
%  plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
end
figure('name','Mode','NumberTitle','off');
imagesc((plot_var))
colorbar
shading interp
set(gcf,'color','w');
axis off
%% Reconstruct velocity field with reduced order
num_modes = 5;
[u_recon] = f_post_POD_reconstruct(var_pod,Sigma,Psy,Phi,num_modes);
%% Show result
figure('name','Original vs reduced order reconstruction','NumberTitle','off');
set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.5]);
for tt=1:options.n_frames
    
plot_var_ori = var_pod(:,:,tt);
plot_var_recon = u_recon(:,:,tt);
if smooth_flag == true
 plot_var = nanconv(plot_var,imageFilter, 'nanout');
end
if exist('im_roi','var')
%  plot_var(im_roi(grid_win:grid_win:size(im_roi,1),grid_win:grid_win:size(im_roi,2)) == 0) = 0;
end
subplot(1,2,1) % original data
imagesc((plot_var_ori))
colorbar
shading interp
caxis([0 2.5])
axis off
subplot(1,2,2) % reconstructed data
imagesc((plot_var_recon))
colorbar
shading interp
caxis([0 2.5])
axis off
set(gcf,'color','w');

drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
