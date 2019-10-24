function [traj_filt,sig_filt] = f_track_filter_trajectories(traj,max_sigma_total,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter of trajectories based on doubled particles and universal outlier
% detection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      traj            - trajectories [Cell Array]
%   ------
%               max_sigma_total - maximum value of tracking criterion
% 
%               options         - options structure
%
%   Output:     traj_filt       - filtered trajectories [Cell Array]
%   -------
%               sig_filt        - tracking criterion of filtered tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 16.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create Filter Trajectory Cell Array
traj_filt = traj;
sig_filt = max_sigma_total;

% Delete Trajectories with same particles
% at t+1
traj_mat = cell2mat(traj);
sig_filt = max_sigma_total;

traj_mat=traj_mat(1:4:end,2);
[n, bin] = histc(traj_mat, unique(traj_mat));
multiple = find(n > 1);
index = find(ismember(bin, multiple));
traj_filt(index) = [];
sig_filt(index) = [];

% at t+2
traj_mat = cell2mat(traj_filt);
traj_mat=traj_mat(1:4:end,3);
[n, bin] = histc(traj_mat, unique(traj_mat));
multiple = find(n > 1);
index = find(ismember(bin, multiple));
traj_filt(index) = [];
sig_filt(index) = [];

% at t+3
traj_mat = cell2mat(traj_filt);
traj_mat=traj_mat(1:4:end,4);
[n, bin] = histc(traj_mat, unique(traj_mat));
multiple = find(n > 1);
index = find(ismember(bin, multiple));
traj_filt(index) = [];
sig_filt(index) = [];

if size(traj_filt,1) > options.n_outlier
    traj_mat = cell2mat(traj_filt);
    x_px=traj_mat(1:4:end,1);
    y_px=traj_mat(2:4:end,1);
    u_px=traj_mat(3:4:end,1);
    v_px=traj_mat(4:4:end,1);

    valid_flag = f_post_outlier_median(x_px,y_px,u_px,v_px,options.n_outlier,options.thr,options.noise);
    n_out = sum(~valid_flag);
    disp([num2str(sum(~valid_flag)),' outliers detected.']);

    traj_filt = traj_filt(valid_flag);
    sig_filt(~valid_flag) = [];

    if size(traj_filt,1) > options.n_outlier
        traj_mat = cell2mat(traj_filt);
        x_px=traj_mat(1:4:end,2);
        y_px=traj_mat(2:4:end,2);
        u_px=traj_mat(3:4:end,2);
        v_px=traj_mat(4:4:end,2);

        valid_flag = f_post_outlier_median(x_px,y_px,u_px,v_px,options.n_outlier,options.thr,options.noise);
        n_out = sum(~valid_flag);
        disp([num2str(sum(~valid_flag)),' outliers detected.']);

        traj_filt = traj_filt(valid_flag);
        sig_filt(~valid_flag) = [];

        if size(traj_filt,1) > options.n_outlier
            traj_mat = cell2mat(traj_filt);
            x_px=traj_mat(1:4:end,3);
            y_px=traj_mat(2:4:end,3);
            u_px=traj_mat(3:4:end,3);
            v_px=traj_mat(4:4:end,3);

            valid_flag = f_post_outlier_median(x_px,y_px,u_px,v_px,options.n_outlier,options.thr,options.noise);
            n_out = sum(~valid_flag);
            disp([num2str(sum(~valid_flag)),' outliers detected.']);

            traj_filt = traj_filt(valid_flag);
            sig_filt(~valid_flag) = [];
        end
    end
end
end

