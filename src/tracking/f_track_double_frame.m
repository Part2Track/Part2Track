function [results, results_px] = f_track_double_frame(part_A_all,part_B_all,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle Tracking of double frame images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      part_A_all       - Particles in image A [n x 5]
%   ------
%               part_B_all       - Particles in image B [n x 5]
%
%               options          - options structure
% 
%
%   Output:     results          - Velocity field [n x 5]
%   -------
%               results_px       - Velocity field in px/time step [n x 5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 16.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image parameters
acq_mode = options.acq_mode; % image acquisition mode
im_res = options.im_res; % camera resolution
n_frames = options.n_frames; % number of double-frames 
dt = options.dt; % time separation of image pairs [s]
m = options.m; % mapping scale [mm/px]
im_roi = options.im_roi; % import image mask, if defined

% particle detection
p_size = options.p_size; %Particle size in px
p_int  = options.p_int; %Particle intensity

% tracking parameter
track_method = options.track_method; % algorithm for particle matching
f_o_s = options.f_o_s; %Field of search in px
n_neighbours = options.n_neighbours; %Number of neigbours for displacement detection
gauss_interp = options.gauss_interp; % 1-gauss fit for hist match, 0-no gauss fit

% outlier detection
n_outlier = options.n_outlier; %Number of neighbours for outlier detection
thr = options.thr; %Threshold for outlier detection
noise = options.noise; %Estimated noise level of measurement

% multi-pass
n_outl_iter = options.n_outl_iter; %Number of iteration for outlier convergence
n_mp = options.n_mp; %Number of steps to increase field of search

% plotting options
plot_int_results = options.plot_int_results; % Show intermediate results: 1-yes 0-no


part_A = part_A_all;
part_B = part_B_all;
%% Main Loop
   
    %% Tracking
    disp('###Perform particle tracking###')
    f_o_s_loop = f_o_s/n_mp;
    k_iter_total = 1;
    while f_o_s_loop <= f_o_s
        k_iter = 1;
        n_out = 0;
        n_matches =1;
        while n_out < n_matches && k_iter <= n_outl_iter
            disp(['Particles left in image A: ',num2str(size(part_A,1))])
            disp(['Particles left in image B: ',num2str(size(part_B,1))])

            switch track_method
                case 'hist_match' 
                   matches = f_track_hist_match( part_A, part_B, f_o_s_loop, n_neighbours, gauss_interp );
                case 'nearest'
                    matches = f_track_nearest_neighbour( part_A, part_B, f_o_s );
            end

            % If matches are found, add to list
            if ~isempty(matches)
                % delete multiple matches
                [n, bin] = histc(matches(:,2), unique(matches(:,2)));
                multiple = find(n > 1);
                index = find(ismember(bin, multiple));
                matches(index,:) = [];
                n_matches = size(matches,1);
                disp(['Number of tracks found: ', num2str(size(matches,1))]);

                if k_iter_total==1
                    x_px = part_A(matches(:,1),1); %[px]
                    y_px = part_A(matches(:,1),2); %[px]
                    x_pxB = part_B(matches(:,2),1); %[px]
                    y_pxB = part_B(matches(:,2),2); %[px]
                    u_px = part_B(matches(:,2),1)-part_A(matches(:,1),1); %[px/dt]
                    v_px = part_B(matches(:,2),2)-part_A(matches(:,1),2); %[px/dt]
                else
                    x_px_temp = part_A(matches(:,1),1); %[px]
                    y_px_temp = part_A(matches(:,1),2); %[px]
                    x_pxB_temp = part_B(matches(:,2),1); %[px]
                    y_pxB_temp = part_B(matches(:,2),2); %[px]
                    u_px_temp = part_B(matches(:,2),1)-part_A(matches(:,1),1); %[px/dt]
                    v_px_temp = part_B(matches(:,2),2)-part_A(matches(:,1),2); %[px/dt]
                    if exist('x_px','var')
                        x_px =[x_px_temp; x_px];
                        y_px =[y_px_temp; y_px];
                        x_pxB =[x_pxB_temp; x_pxB];
                        y_pxB =[y_pxB_temp; y_pxB];
                        u_px =[u_px_temp; u_px];
                        v_px =[v_px_temp; v_px];
                    else
                        x_px = x_px_temp;
                        y_px = y_px_temp;
                        x_pxB = x_pxB_temp;
                        y_pxB = y_pxB_temp;
                        u_px = u_px_temp;
                        v_px = v_px_temp;
                    end
                end

            %% Outlier Detection    
            disp('###Perform outlier detection###')

            valid_flag = f_post_outlier_median(x_px,y_px,u_px,v_px,n_outlier,thr,noise);
            n_out = sum(~valid_flag);
            disp([num2str(sum(~valid_flag)),' outliers detected.']);

            x_px_inv = x_px(~valid_flag);
            x_pxB_inv = x_pxB(~valid_flag);
            x_px = x_px(valid_flag);
            y_px = y_px(valid_flag);
            x_pxB = x_pxB(valid_flag);
            y_pxB = y_pxB(valid_flag);
            u_px = u_px(valid_flag);
            v_px = v_px(valid_flag);
            %% Calculate Coordinates and Velocities
            x = x_px*m; %[mm]
            y = y_px*m; %[mm]
            u = (x_pxB-x_px)*m/dt/1000; %[m/s]
            v = (y_pxB-y_px)*m/dt/1000; %[m/s]
            valid_flag(valid_flag==0)=[];
            
            % save matches to result matrix
            results_px = [x_px y_px u_px v_px valid_flag];
            results = [x y u v valid_flag];

            % plot intermediate results
            if plot_int_results == 1
                if k_iter_total==1
                    g=figure('name','Intermediate Results','NumberTitle','off','Color','w');
                    set(g,'Units', 'centimeter', 'Position', [5, 5, 10, 10])
                else
                    if exist('g','var')
                    clf(g)
                    end
                end
                scatter(part_A_all(:,1)*m,part_A_all(:,2)*m,2,'filled')
                hold on
                quiver(x(valid_flag),y(valid_flag),u(valid_flag),v(valid_flag),1.0,'black','linewidth',1.0)    
                xlim([0 im_res(2)*m])
                ylim([0 im_res(1)*m])
                axis off
                axis equal
                set(gca,'ydir','reverse')    
                pause(0.1)
            end 

            [ind_log,~] = ismember(part_A_all(:,1),x_px);    
            part_A = part_A_all;
            part_A(ind_log,:) =[];
            [~,ind_uni] = unique(part_A(:,1));
            part_A = part_A(ind_uni,:);                      
            [ind_log,~] = ismember(part_B_all(:,1),x_pxB);    
            part_B = part_B_all;
            part_B(ind_log,:) =[];
            [~,ind_uni] = unique(part_B(:,1));
            part_B = part_B(ind_uni,:);
            
            end

            clear x y u v matches

            if n_out >= n_matches || k_iter == n_outl_iter
                f_o_s_loop = f_o_s_loop + f_o_s/n_mp;
            end
            k_iter = k_iter+1;
            k_iter_total = k_iter_total+1;
        end
    end
end

