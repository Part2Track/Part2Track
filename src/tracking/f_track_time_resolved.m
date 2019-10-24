function [traj] = f_track_time_resolved(dir_save,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-resolved processing of image sequence with four-frame trajectory
% initialization and extrapolation in next time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      dir_save      - Path to save folder [string]
%   ------
%               options       - options structure
% 
%
%   Output:     traj          - trajectories [Cell Array]
%   -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 16.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_dist = options.min_dist;

%% Initialize first particle trajectories
part_all{4,1} = [];
for tt=1:4 %
    load([dir_save,'part_',num2str(tt),'.mat'])
    part_all{tt}=part;
end

%% Create global velocity field for first time step
[~, results_px] = f_track_double_frame(part_all{1},part_all{1+1},options);
x = results_px(:,1);
y = results_px(:,2);
u = results_px(:,3);
v = results_px(:,4);    

Fu = scatteredInterpolant(x,y,u);
Fv = scatteredInterpolant(x,y,v);

%% Trajectory reconstruction for first four time steps
part_res = part_all;
for kk=1:options.n_mp_ti % Iterative Trajectory Reconstruction
    disp('###Initialize Trajectories###')
    [traj,max_sigma_total] = f_track_init_trajectories(part_res,min_dist,Fu,Fv,options);
    
    disp(['Trajectories found: ',num2str(size(traj,1))])
    if kk==1
        traj_final = traj;
        max_sigma_final = max_sigma_total;
    else
        traj_final = [traj_filt; traj];
        max_sigma_final = [sig_filt max_sigma_total];
    end

    % Filter trajectories for outliers
    disp('###Filter Trajectories###')
    [traj_filt,sig_filt] = f_track_filter_trajectories(traj_final,max_sigma_final,options);

    %% Delete correct particles from particle list
    part_res = part_all;
    for jj = 1:4
        traj_filt_mat = cell2mat(traj_filt);
        part_temp = [traj_filt_mat(1:4:end,jj) traj_filt_mat(2:4:end,jj)];

        index_mem = ismember(part_res{jj}(:,1:2), part_temp);
        index_mem = (index_mem(:,1) == 1 & index_mem(:,2) == 1);

        part_res{jj}(index_mem,:) = [];
        
    end
end
traj = traj_filt;
max_sigma_final = sig_filt;
% Plot results
if options.plot_int_results == 1
    [~,h] = f_post_plot_trajectories(traj,options);
    pause(0.1)
end

%% Advance Trajectories in time
for tt = 5:options.n_frames
    load([dir_save,'part_',num2str(tt),'.mat'])
    part_all{tt}=part;
    part_res{tt}=part;

    disp(['###Timestep: ',num2str(tt),'###'])
    disp('###Extrapolate Trajectories###')
    % Loop over all tracked particles
    for i_part = 1:size(traj,1)
        if isnan(traj{i_part}(1,tt-1)) == 0
            xfit = polyfit([1:4]',traj{i_part}(1,tt-4:tt-1)',2);
            yfit = polyfit([1:4]',traj{i_part}(2,tt-4:tt-1)',2);
            x_pred = polyval(xfit,5);
            y_pred = polyval(yfit,5);

            dist = sqrt((part_res{tt}(:,1)-x_pred).^2+(part_res{tt}(:,2)-y_pred).^2);
            [min_dist_pred,i_dist] = min(dist);
            min_dist_save(i_part) = min_dist_pred;
            
            if abs(min_dist_pred) < min_dist
                traj{i_part}(1,tt) = part_res{tt}(i_dist,1);
                traj{i_part}(2,tt) = part_res{tt}(i_dist,2);
                traj{i_part}(3,tt) = part_res{tt}(i_dist,1)-traj{i_part}(1,tt-1);
                traj{i_part}(4,tt) = part_res{tt}(i_dist,2)-traj{i_part}(2,tt-1);
                part_res{tt}(i_dist,:) = [];
            else
                traj{i_part}(1,tt) = NaN;
                traj{i_part}(2,tt) = NaN;
                traj{i_part}(3,tt) = NaN;
                traj{i_part}(4,tt) = NaN;
            end
        else
            traj{i_part}(1,tt) = NaN;
            traj{i_part}(2,tt) = NaN;
            traj{i_part}(3,tt) = NaN;
            traj{i_part}(4,tt) = NaN;
        end            
    end

    % Find new trajectories among residual particles
    disp('###Reconstruct new Trajectories###')
    [traj_new,max_sigma_total] = f_track_init_trajectories(part_res(tt-3:tt),min_dist,[],[],options);
    disp(['Trajectories found: ',num2str(size(traj_new,1))])
    
    disp('###Filter new Trajectories###')
    if ~isempty(traj_new) == 1
        [traj_new,sig_filt] = f_track_filter_trajectories(traj_new,max_sigma_total,options);

        % Add NaNs to previous not tracked timesteps
        traj_mat = cell2mat(traj_new);
        traj_mat = [NaN(size(traj_mat,1),tt-size(traj_mat,2)) traj_mat];
        traj_new = mat2cell(traj_mat,[4*ones(size(traj_new,1),1)], [tt]);

        % Add new trajectories to existing list
        traj = [traj; traj_new];
        max_sigma_final = [max_sigma_final sig_filt];
        %% Filter short trajectories

        traj_mat = cell2mat(traj);
        index_del = (sum(~isnan(traj_mat')) < 5) & (isnan(traj_mat(:,tt)'));
        traj_mat(index_del,:) = [];        
        max_sigma_final(index_del(1:4:end)) = [];        
        traj = mat2cell(traj_mat,[4*ones(size(traj_mat,1)/4,1)], [size(traj_mat,2)]);
        
        assignin('base','traj',traj)

        % Update Plot
        if options.plot_int_results == 1 
            traj_mat = cell2mat(traj);
            X_tr_plot = traj_mat(1:4:end,:)';
            Y_tr_plot = traj_mat(2:4:end,:)';
            disp_plot = sqrt(traj_mat(3:4:end,:).^2+traj_mat(4:4:end,:).^2)';
            set(h,'XData',[X_tr_plot ; nan(1,size(X_tr_plot,2))],'YData',[Y_tr_plot ; nan(1,size(Y_tr_plot,2))],'CData',[(disp_plot); nan(1,size(X_tr_plot,2))]);
            pause(0.1)
        end
    end
end

