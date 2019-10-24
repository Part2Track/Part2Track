function [u_recon] = f_post_POD_reconstruct(u_pod,Sigma,Psy,Phi,num_modes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruct given field u_pod with a defined number of POD modes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      u_pod     - previous decomposed variable
%   ------
%               Psy       - Left singular vectors / Spatial modes
%   
%               Sigma     - Singular values / Energy content
%
%               Phi       - Right singular vector / Temporal information
%
%               num_modes - number of modes used for reconstruction
%
%   Output:     u_recon   - reconstructed field     
%   -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 26.07.19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Reconstruct velocity field with defined order
for i_mode = 1:num_modes
    if i_mode == 1
        u_recon = Sigma(i_mode)*Psy(:,i_mode)*Phi(:,i_mode)';
    else
        u_recon = u_recon + Sigma(i_mode)*Psy(:,i_mode)*Phi(:,i_mode)';
    end
end
u_recon = reshape(u_recon,size(u_pod,1),size(u_pod,2),size(u_pod,3));

end

