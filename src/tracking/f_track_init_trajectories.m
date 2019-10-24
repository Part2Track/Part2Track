function [traj,max_sigma_total] = f_track_init_trajectories(part_all,min_dist,Fu,Fv,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of trajectories from particle positions of four time steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      part_all        - Particle coordinates [Cell Array 4x1]
%   ------
%               min_dist        - allowed distance from predicted position
% 
%               Fu              - scattered interpolant of velocity field u
%
%               Fv              - scattered interpolant of velocity field v
%
%               options         - Global options structure
%
%   Output:     traj            - trajectories [Cell Array]
%   -------
%               max_sigma_total - maximum value of tracking criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 16.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_o_s = options.f_o_s;

% Initialize Variables
traj{size(part_all{1}(:,1),1),1} = [];
max_sigma_total=zeros(1,1);

% Loop over all particles
parfor i_part = 1:size(part_all{1},1)
   [traj_temp, sigma_temp]  = f_track_init_trajectories_loop(i_part,part_all,min_dist,Fu,Fv,f_o_s);
   
   % Save found trajectory in final variables
   traj{i_part} = traj_temp;
   max_sigma_total(i_part) = sigma_temp;
end

% Delete empty entries
traj = traj(~cellfun('isempty',traj));
max_sigma_total(isnan(max_sigma_total)) = [];
% max_sigma_total = max_sigma_total';

end
