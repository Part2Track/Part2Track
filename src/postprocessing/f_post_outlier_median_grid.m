function [ valid_flag ] = f_post_outlier_median_grid(u,v,b,thr,noise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universal outlier detection method based on:
% J. Westerweel, F. Scarano: "Univsersal outlier detection for PIV data"
% Experiments in Fluids, 2005, 39:1096 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      x             - x-coordinate of vectors [n x 1]
%   ------
%               y             - y-coordinate of vectors [n x 1]
% 
%               u             - x-velocity-component of vector [n x 1]
% 
%               v             - y-velocity-component of vector [n x 1]
%
%               b             - Kernel-width [integer]
%
%               thr           - fluctuation threshold for 
%                               outlier detection [integer]
%                               |
%                               ---> recommended: thr = 2
%
%               noise          - noise level of measurement [integer]
%                                |
%                                ---> recommended: noise = 0.1
%
%
%   Output:     valid_flag     - logical vector of outlier detection choice
%   -------                         |
%                                   |--> 1: valid vector
%                                   |
%                                   |--> 0: invalid vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 17.11.17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[J,I] = size(u);

median_res = zeros(J,I);
norm_fluc = zeros(J,I);

for cc=1:2
    if cc==1; vel_c=u; else; vel_c=v; end
    for ii=1+b:I-b
        for jj=1+b:J-b
            neigh = vel_c(jj-b:jj+b,ii-b:ii+b);
            neighcol = neigh(:);
            neighcol2 = [neighcol(1:(2*b+1)*b+b); neighcol((2*b+1)*b+b+2:end)];
            
            med=median(neighcol2);
            fluc = vel_c(jj,ii)-med;
            res = neighcol2-med;
            median_res = median(abs(res));
            norm_fluc(jj,ii,cc) =abs(fluc/(median_res+noise));
        end
    end
end
valid_flag = sqrt(norm_fluc(:,:,1).^2 + norm_fluc(:,:,2).^2) < thr; %detection criterion

end

