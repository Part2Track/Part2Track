function [Psy, Sigma, Phi] = f_post_POD_decompose(u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POD decomposition of input field u using singular value decomposition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      u         - field to decompose [ m x n x t]
%   ------
%
%   Output:     Psy       - Left singular vectors / Spatial modes
%   -------     
%               Sigma     - Singular values / Energy content
%
%               Phi       - Right singular vector / Temporal information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 26.07.19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Construct strucute matrix [m*n x t]
u_pod = reshape(u,size(u,1)*size(u,2),size(u,3));

% Set invalid values to zero
u_pod(isnan(u_pod)) = 0;
u_pod(isinf(u_pod)) = 0;

%% Perform singular value decompisition
[Psy,Sigma,Phi] = svd(u_pod,0);
Sigma = diag(Sigma);
end

