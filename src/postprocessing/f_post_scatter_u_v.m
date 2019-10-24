function [fig_handle,scat_handle] = f_post_scatter_u_v(plot_u,plot_v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scatter plot of velocity components u and v color-coded with magnitude of 
% velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      plot_u      - velocity component u
%   ------
%               plot_v      - velocity component v
%
%   Output:     fig_hande   - handle to figure
%   -------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% T. Janke / 26.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis_font = 'Arial'; % Font
axis_font_size = 16; % Fontsize

fig_handle=figure('name','Scatter Plot Velocity','NumberTitle','off');
set(fig_handle, 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
set(gcf,'color','w'); % set background color 

scat_handle=scatter(plot_u(:),plot_v(:),8,sqrt(plot_u(:).^2+plot_v(:).^2),'filled');
xlabel('u');
ylabel('v');
ax1 = gca;
ax1.XAxisLocation = 'origin';
ax1.YAxisLocation = 'origin';
ax1.FontSize=axis_font_size;
ax1.FontName = axis_font;
end

