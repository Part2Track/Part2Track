function [traj,max_sigma_total] = f_track_init_trajectories_loop(i_part,part_all,min_dist,Fu,Fv,f_o_s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of trajectories from particle positions of four time steps
% - Core Function to loop over - 
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
%               f_o_s           - Field of Search
%
%   Output:     traj            - trajectories [Cell Array]
%   -------
%               index_save      - index of saved matches
%
%               max_sigma_total - maximum value of tracking criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 16.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(Fu) == 1 % If no predictor is given, use nearest neighbour
    dist = sqrt((part_all{1+1}(:,1)-part_all{1}(i_part,1)).^2+(part_all{1+1}(:,2)-part_all{1}(i_part,2)).^2);
    indice{1} = find(dist < f_o_s);
else % Use predictor   
    x_pred = part_all{1}(i_part,1) + Fu(part_all{1}(i_part,1),part_all{1}(i_part,2));
    y_pred = part_all{1}(i_part,2) + Fv(part_all{1}(i_part,1),part_all{1}(i_part,2));

    dist = sqrt((part_all{1+1}(:,1)-x_pred).^2+(part_all{1+1}(:,2)-y_pred).^2);
    indice{1} = find(dist < min_dist);
end

if ~isempty(indice{1}) 
    
    % Predict new position at t+1
	coord_pred_1 = (2*part_all{1+1}(indice{1}(:),1:2) - part_all{1}(i_part,1:2))';
    % Find particles with are within min_dist to predicted position
    dist = sqrt(sum(part_all{1+2}(:,1:2).^2,2)+sum(coord_pred_1.^2)-2*part_all{1+2}(:,1:2)*coord_pred_1);
    indice{2} = dist < min_dist;

    % Loop over canditates at t+1
    for i_cand1 = 1:size(indice{1},1)

        if max(indice{2}(:,i_cand1))>0 % Catch if there are no candidates

            % Candidate list
            part_cand =  part_all{1+2}(indice{2}(:,i_cand1),1:2);

            % Build matrices for linear regression
            % X_reg - time steps 1:3
            X_reg = repmat([ones(length([1:3]'),1) [1:3]' [ones(length([1:3]'),1) [1:3]']], 1,size(part_cand,1));

            % y_reg - x- and y-coordinates at t=1:3
            y_reg = [part_all{1}(i_part,1) part_all{1+1}(indice{1}(i_cand1),1);...
                    part_all{1}(i_part,2) part_all{1+1}(indice{1}(i_cand1),2)];
            y_reg = repmat(y_reg, size(part_cand,1),1);
            y_append = part_cand';
            y_append=y_append(:);
            y_reg = [y_reg y_append];

            % Perform linear regression
            coeff_fit = X_reg\y_reg';

            % Predict new position at t+2
            coord_pred_temp = 4*coeff_fit(2,:)+coeff_fit(1,:);
            coord_pred = [coord_pred_temp(1:2:end)' coord_pred_temp(2:2:end)']';

            % Find particles with are within min_dist to predicted position
            dist = sqrt(sum(part_all{1+3}(:,1:2).^2,2)+sum(coord_pred.^2)-2*part_all{1+3}(:,1:2)*coord_pred);
            dist = dist < min_dist;

            % Connect all candidates to possible trajectories
            if max(dist(:)) > 0  % Catch if there are no candidates 
                % Loop over canditates at t+2
                for i_cand2 = 1:size(part_cand,1)
                    num_cand = sum(dist(:,i_cand2));
                    traj_append = zeros(num_cand*2,4);
                    traj_append(1:2:end) = [repmat([part_all{1}(i_part,1) part_all{1+1}(indice{1}(i_cand1),1) part_cand(i_cand2,1)],num_cand,1) part_all{1+3}(dist(:,i_cand2),1)];
                    traj_append(2:2:end) = [repmat([part_all{1}(i_part,2) part_all{1+1}(indice{1}(i_cand1),2) part_cand(i_cand2,2)],num_cand,1) part_all{1+3}(dist(:,i_cand2),2)];

                    if ~exist('traj_temp','var')
                        traj_temp = traj_append;
                    else
                        traj_temp = [traj_temp; traj_append];
                    end

                end
            end
        end
    end
end
% traj = cell2mat(traj_temp);
% Find best trajectory under all possible trajectories
if exist('traj_temp','var') 
    kk=size(traj_temp,1)+1;
    velo = diff(traj_temp')';
    L = sqrt(velo(1:2:end,:).^2+velo(2:2:end,:).^2);
    theta = zeros((kk-1)/2,2);
    sigma_L = zeros((kk-1)/2,1);
    sigma_theta = zeros((kk-1)/2,1);
    sigma_total = zeros((kk-1)/2,1);
    for jj=1:(kk-1)/2
        d1 = [traj_temp(jj*2-1,2)-traj_temp(jj*2-1,1);traj_temp(jj*2,2)-traj_temp(jj*2,1); 0];
        d2 = [traj_temp(jj*2-1,3)-traj_temp(jj*2-1,2);traj_temp(jj*2,3)-traj_temp(jj*2,2); 0];
        theta(jj,1) = atan2d(norm(cross(d1,d2)),dot(d1,d2));

        d1 = [traj_temp(jj*2-1,3)-traj_temp(jj*2-1,2);traj_temp(jj*2,3)-traj_temp(jj*2,2); 0];
        d2 = [traj_temp(jj*2-1,4)-traj_temp(jj*2-1,3);traj_temp(jj*2,4)-traj_temp(jj*2,3); 0];
        theta(jj,2) = atan2d(norm(cross(d1,d2)),dot(d1,d2));

        sigma_L(jj,1)=sqrt(((L(jj,1)-mean(L(jj,:)))^2+(L(jj,2)-mean(L(jj,:)))^2+(L(jj,3)-mean(L(jj,:)))^2)/3);
        sigma_theta(jj,1)=sqrt(((theta(jj,1)-mean(theta(jj,:)))^2+(theta(jj,2)-mean(theta(jj,:)))^2+(L(jj,3)-mean(L(jj,:)))^2)/3);

        sigma_total(jj,1) = sqrt(sigma_L(jj,1).^2/mean(L(jj,:)).^2+sigma_theta(jj,1).^2);
    end

    final_indice = find( sigma_total == min(sigma_total));

    % in the case of no particle motion
    if isempty(final_indice) == 1 && isnan(min(sigma_total)) == 1
        final_indice = 1;
        sigma_total = 99999;
    end

    traj= [traj_temp(final_indice(1),:);traj_temp(final_indice(1)+1,:);...
                    [velo(final_indice(1),:) velo(final_indice(1),end)];...
                    [velo(final_indice(1)+1,:) velo(final_indice(1)+1,end)]];
    max_sigma_total = max(sigma_total(final_indice(1)));

else
    traj = [];
    max_sigma_total = NaN;
end

%% Non vectorized code - just in case ;)

% kk = 1;    
% 
% if isempty(Fu) == 1 % If no predictor is given, use nearest neighbour
%     dist = sqrt((part_all{1+1}(:,1)-part_all{1}(i_part,1)).^2+(part_all{1+1}(:,2)-part_all{1}(i_part,2)).^2);
%     indice{1} = find(dist < f_o_s);
% else % Use predictor   
%     x_pred = part_all{1}(i_part,1) + Fu(part_all{1}(i_part,1),part_all{1}(i_part,2));
%     y_pred = part_all{1}(i_part,2) + Fv(part_all{1}(i_part,1),part_all{1}(i_part,2));
% 
%     dist = sqrt((part_all{1+1}(:,1)-x_pred).^2+(part_all{1+1}(:,2)-y_pred).^2);
%     indice{1} = find(dist < min_dist);
% end
% 
% % Loop over canditates at t+1
% for i_cand1 = 1:size(indice{1},1)
% 
%     x_pred = 2*part_all{1+1}(indice{1}(i_cand1),1) - part_all{1}(i_part,1);
%     y_pred = 2*part_all{1+1}(indice{1}(i_cand1),2) - part_all{1}(i_part,2);
% 
%     dist = sqrt((part_all{1+2}(:,1)-x_pred).^2+(part_all{1+2}(:,2)-y_pred).^2);
%     indice{2} = find(dist < min_dist);
% 
%     % Loop over canditates at t+2
%     for i_cand2 = 1:size(indice{2},1)
% 
% 
%         x_fit = polyfit([1:3]',[part_all{1}(i_part,1) part_all{1+1}(indice{1}(i_cand1),1) part_all{1+2}(indice{2}(i_cand2),1)]',1);
%         y_fit = polyfit([1:3]',[part_all{1}(i_part,2) part_all{1+1}(indice{1}(i_cand1),2) part_all{1+2}(indice{2}(i_cand2),2)]',1);
%    
%         x_pred = 4*x_fit(1)+x_fit(2);
%         y_pred = 4*y_fit(1)+y_fit(2);
% 
%         dist = sqrt((part_all{1+3}(:,1)-x_pred).^2+(part_all{1+3}(:,2)-y_pred).^2);
%         indice{3} = find(dist < min_dist);
%         
%         % Loop over canditates at t+3
%         for i_cand3 = 1:size(indice{3},1)
%             traj_temp(kk,1:4) = [part_all{1}(i_part,1) part_all{1+1}(indice{1}(i_cand1),1) part_all{1+2}(indice{2}(i_cand2),1) part_all{1+3}(indice{3}(i_cand3),1)];
%             traj_temp(kk+1,1:4) = [part_all{1}(i_part,2) part_all{1+1}(indice{1}(i_cand1),2) part_all{1+2}(indice{2}(i_cand2),2) part_all{1+3}(indice{3}(i_cand3),2)];
%             index_temp(kk,1:4) =[i_part, indice{1}(i_cand1) indice{2}(i_cand2) indice{3}(i_cand3)];
%             kk=kk+2;
%         end
%     end
% end
end

