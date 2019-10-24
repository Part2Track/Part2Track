function [Delta,Lambda_ci,Q,Lambda_2] = f_post_vortex_criteria(u,v,dh,scale_rms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of four different vortex criteria. The calculation is based
% on:
% Q. Chen, Q. Zhong, M. Qi, X. Wang: "Comparison of vortex identification
% criteria for planaer velocity fields in wall turbulence"
% Phys. Fluids 27, 085101(2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      u         - x - velocity component
%   ------
%               v         - y - velocity component 
% 
%               dh        - uniform vector spacing
%
%               scale_rms - scale with rms ('true' or 'false')
%
%   Output:     Delta     - Delta Criterion
%   -------     
%               Lambda_ci - Lambda_ci Criterion
%
%               Q         - Q Criterion
%
%               Lambda_2  - Lambda_2 Criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 10.05.19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Input check
if nargin < 2
  error('u & v velocity fields missing!')
end

if nargin < 3
  dh = 1;
end

if nargin < 4
  scale_rms = 'false';
end

% Calculate spatial gradients
[dudx,dudy] = gradient(u,dh);
[dvdx,dvdy] = gradient(v,dh);

% Delta Criterion
% Connected region of negative Delta -> Vortex Identification
Delta = 4*dudy.*dvdx + (dudx - dvdy).^2;

% Lambda_ci Criterion
Lambda_ci = 1/2*sqrt( -4*(dudy.*dvdx-dudx.*dvdy) - (dudx+dvdy).^2);
% Connected regions of non-zero imaginary part of Lambda_ci -> Vortex
% Identification
Lambda_ci = imag(Lambda_ci);

% Q Criterion
% Connected regions of positive Q -> Vortex Identification
Q = -dudy.*dvdx - 1/2*dudx.^2 - 1/2*dvdy.^2;

% Lambda_2 Criterion
% Connected regions of negative Lambda_2 -> Vortex Identification
Lambda_2 = dudy.*dvdx + 1/2*( dudx.^2+dvdy.^2 + abs(dudx+dvdy).*sqrt( (dudx-dvdy).^2 + (dudy+dvdx).^2));

% Scale with rms(Criterion) if selected
if scale_rms == true
    Delta       = Delta./rms(rms(Delta));
    Lambda_ci	= Lambda_ci./rms(rms(Lambda_ci));
    Q           = Q./rms(rms(Q));
    Lambda_2	= Lambda_2./rms(rms(Lambda_2));
end
end

