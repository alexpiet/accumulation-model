function [] = plot_pdf(D)
% Plots the PDF at each time point
% D.plot_zero_line, will plot a dashed line at a=0
% D.plot_mean_line, will plot a line at the mean of the PDF at each timepoint

if ~isfield(D, 'zlimit'); D.zlimit = 0.5; end;
imagesc(D.avals, D.T, D.pdf, [0 .5])
colormap hot;
ylabel('Time (s)', 'fontsize',12)
xlabel('Accumulated Evidence (a)' ,'fontsize',12)
set(gca,'fontsize',12);

if ~isfield(D, 'plot_zero_line'); D.plot_zero_line =1; end;
if D.plot_zero_line
    hold on;
    plot([0 0], [D.T(1) D.T(end)],'m--','linewidth',2)
end

if ~isfield(D, 'plot_mean_line'); D.plot_mean_line =0; end;
if D.plot_mean_line
    hold on;
    if isfield(D, 'ma')
        plot(D.ma, D.T,'b-','linewidth',2)
    elseif isfield(D, 'mean')
        plot(D.mean, D.T,'b-','linewidth',2)
    end
end

if isfield(D, 'left_clicks');
    plot(-D.click_lim, D.left_clicks, '>','color',[48 127 255]./255,'markerfacecolor',[48 127 255]./255,'markersize',D.click_size)
end
if isfield(D, 'right_clicks');
    plot(D.click_lim, D.right_clicks, '<','color',[0 140 54]./255 ,'markerfacecolor',[0 140 54]./255,'markersize',D.click_size)
end

if isfield(D, 'title')
    title(D.title)
end
xlim([D.avals(1) D.avals(end)])
ylim([D.T(1) D.T(end)])



