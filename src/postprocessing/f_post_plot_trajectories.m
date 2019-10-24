function [f,h] = f_post_plot_trajectories(traj,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots all trajectories color-coded with the magnitude of velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      traj            - trajectories [Cell Array]
%   ------
%               options         - options structure
%
%   Output:     f               - figure handle
%   -------
%               h               - patch handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 16.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im_res = options.im_res;
m = options.m;

%% Plot intermediate results
traj_mat = cell2mat(traj);
X_tr_plot = traj_mat(1:4:end,:)';
Y_tr_plot = traj_mat(2:4:end,:)';
disp_plot = sqrt(traj_mat(3:4:end,:).^2+traj_mat(4:4:end,:).^2)';
f=figure('name','Trajectories','NumberTitle','off','Color','w');

h=patch([X_tr_plot ; nan(1,size(X_tr_plot,2))],[Y_tr_plot; nan(1,size(X_tr_plot,2))],[(disp_plot); nan(1,size(X_tr_plot,2))],'EdgeColor','interp','FaceColor','none','LineWidth',1.0);
set(gca,'YDir','reverse')

end

