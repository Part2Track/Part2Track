function [fig_handle] = f_post_check_peak_locking(plot_var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot distribution of particle displacements in order to evaluate possible
% peak locking effects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      plot_var    - variable to plot
%   ------
%
%   Output:     fig_hande   - handle to figure
%   -------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% T. Janke / 26.11.18
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_var(isnan(plot_var))=[]; % Delete NaNs
plot_var(plot_var==0) = []; % Delete Zeros
plot_var_sub = (abs(plot_var)-floor(abs(plot_var)))-0.5; % Calculate subpixel shifts

axis_font = 'Arial'; % Font
axis_font_size = 16; % Fontsize

fig_handle=figure('name','Peak Locking Check','NumberTitle','off');
set(fig_handle, 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);
set(gcf,'color','w'); % set background color    

subplot(2,1,1)
% Histrogram of shift
histogram(plot_var(:),min(plot_var(:)):0.1:max(plot_var(:)));
    set(gca,'FontSize',axis_font_size,'FontName',axis_font);
    box off
    xlabel('Particle Displacement [px]')
    xlim([min(plot_var(:)) max(plot_var(:))])
    ylabel('Counts')
subplot(2,1,2)
% Histogram of sub-pixel shift
h1=histogram(plot_var_sub,-0.5:0.05:0.5);
    set(gca,'FontSize',axis_font_size,'FontName',axis_font);
    box off
    xlabel('Particle Subpixel Displacement [px]')
    xlim([-0.5 0.5])
    ylabel('Counts')
    
% Calculate Peak locking degree
N_min = min(h1.BinCounts);
N_max = max(h1.BinCounts);
C = 1-N_min/N_max;

text(-0.3,max(h1.BinCounts),['Peak locking degree: ',num2str(C)],...
    'FontSize',axis_font_size','FontName',axis_font);
end

