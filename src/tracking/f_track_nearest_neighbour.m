function [ matches ] = f_track_nearest_neighbour( part_A, part_B, f_o_s )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nearest neigbhour tracking routine.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      part_A      - coordinates of particles in image A [n x 4]
%   ------
%               part_B      - coordinates of particles in image B [n x 4]
% 
%               f_o_s       - field of search [px]
%
%   Output:     matches     - list of indices of matching particles [m x 2]
%   -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 21.09.17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk = 1;
for ii=1:size(part_A,1) % Loop all particles in Image A
    
    temp_dist = sqrt((part_B(:,1)-part_A(ii,1)).^2+(part_B(:,2)-part_A(ii,2)).^2); % Calculate all possible distances between particle ii and particles in image B
    
    if min(temp_dist) < f_o_s % Just consider candidates within field of search
        [~,index_min] = min(temp_dist); % Find candidate with minimum displacement
        matches(kk,1:2) = [ii, index_min];
        kk = kk+1;
    end
end

end

