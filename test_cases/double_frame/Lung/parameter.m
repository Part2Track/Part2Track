%% recording parameters
options.acq_mode = 'double_frame'; %parameter.acq_mode; % image acquisition mode
options.im_res = [1200 1600]; %parameter.im_res; % camera resolution
options.n_frames = 10; % parameter.n_frames; % number of double-frames 
options.dt = 1.6e-4; %parameter.dt; % time separation of image pairs [s]
options.m = 0.0408; %parameter.m; % mapping scale [mm/px]
options.im_pre = 'run_'; % image file prefix
options.im_file = 'tif'; % image file format
load('im_roi.mat'); options.im_roi = im_roi; % import image mask, if defined
%options.im_roi = ones(options.im_res(1),options.im_res(2)); % define mask else

%% Processing Parameters
% particle detection
options.p_size = 3; %Particle size in px
options.p_int  = 3000; %Particle intensity

%% Double-Frame processing (also utilized in Time-Resolved Processing)
% tracking parameter
options.track_method = 'hist_match'; % algorithm for particle matching
options.f_o_s = 25; %Field of search in px
options.n_neighbours = 25; %Number of neigbours for displacement detection
options.gauss_interp = 1; % 1-gauss fit for hist match, 0-no gauss fit

% outlier detection
options.n_outlier = 25; %Number of neighbours for outlier detection
options.thr = 2; %Threshold for outlier detection
options.noise = 0.1; %Estimated noise level of measurement

% multi-pass
options.n_outl_iter = 3; %Number of iteration for outlier convergence
options.n_mp = 3; %Number of steps to increase field of search

%% Time-resolved processing
options.min_dist = 1; %Maximum deviation from position prediction in px
options.n_mp_ti = 1; %Number of iteration for first track initialization