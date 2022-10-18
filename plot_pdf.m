function [] = plot_pdf(D)
% Plots the PDF at each time point
% D.plot_zero_line, will plot a dashed line at a=0
% D.plot_mean_line, will plot a line at the mean of the PDF at each timepoint

if ~isfield(D, 'left_color')
    D.left_color = [48 127 255]./255;
end
if ~isfield(D, 'right_color')
    D.right_color = [0 140 54]./255 ;
end
left_color = D.left_color;
right_color = D.right_color;


if ~isfield(D, 'zlimit'); D.zlimit = 0.15; end;
%imagesc(D.avals, D.T, D.pdf, [0 .5])

imagesc(D.T,D.avals,D.pdf',[0 D.zlimit])
hold on
xlabel('Time (s)', 'fontsize',12)
ylabel('Accumulated Evidence (a)' ,'fontsize',12)
%%
if ~isfield(D, 'plot_zero_line'); D.plot_zero_line =1; end;
if D.plot_zero_line
    hold on;
    %plot([0 0], [D.T(1) D.T(end)],'m--','linewidth',2)
    plot([D.T(1) D.T(end)],[0 0],'color',[1 1 1].*.0,'linewidth',1)
end

if ~isfield(D, 'plot_bias_line'); D.plot_bias_line =0; end;
if D.plot_bias_line
    hold on;
    %plot([0 0], [D.T(1) D.T(end)],'m--','linewidth',2)
    plot([D.T(1) D.T(end)],[0 0]+D.bias_param,'--',...
        'color',[1 1 1].*.0,'linewidth',1)
end

if ~isfield(D, 'plot_mean_line'); D.plot_mean_line =0; end;
if D.plot_mean_line
    hold on;
    if isfield(D, 'ma')
        %plot(D.ma, D.T,'b-','linewidth',2)
        plot(D.T, D.ma,'k-','linewidth',2)
    elseif isfield(D, 'mean')
        %plot(D.mean, D.T,'b-','linewidth',2)
        plot(D.T, D.mean,'k-','linewidth',2)
    end
end
%%
if isfield(D,'model_switches')
    plot([D.model_switches; D.model_switches],D.model_switch_y-[0;1],...
        '-','color', 'r','linewidth',2)
end
if isfield(D,'state_switches')
    plot([D.state_switches; D.state_switches],D.state_switch_y-[0;1],...
        '-','color', 'k','linewidth',2)
end
%%
if isfield(D,'left_click_marker') & isempty(D.left_click_marker)
    D.left_click_marker= '<';
end
if isfield(D,'right_click_marker') & isempty(D.right_click_marker)
    D.right_click_marker= '>';
end



if isfield(D, 'left_clicks');
    if ~isfield(D,'left_click_y') & isfield(D,'click_lim')
        D.left_click_y = -D.click_lim;
    end
    if D.left_click_marker == '|'
        nl = length(D.left_clicks);
        xx = [D.left_clicks; D.left_clicks]';
        yy =  D.left_click_y-repmat([0 1],nl,1).*D.click_height;
        plot(xx', yy',  '-', 'color', left_color, 'linewidth', 1.5)
    else
        plot(D.left_clicks, D.left_click_y,  D.left_click_marker, ...
            'color',left_color,'markerfacecolor','w',...
            'markersize',D.click_size)
    end
    
end
if isfield(D, 'right_clicks');
    if ~isfield(D,'right_click_y') & isfield(D,'click_lim')
        D.right_click_y = D.click_lim;
    end
    if D.right_click_marker == '|'
        nr = length(D.right_clicks);
        xx = [D.right_clicks; D.right_clicks]';
        yy =  D.right_click_y+repmat([0 1],nr,1).*D.click_height;
        plot(xx', yy',  '-', 'color', right_color, 'linewidth', 1.5)
    else
        plot(D.right_clicks, D.right_click_y,  D.right_click_marker, ...
            'color',right_color,'markerfacecolor','w',...
            'markersize',D.click_size)
    end
end

if isfield(D, 'title')
    title(D.title)
end
ylim([D.avals(1) D.avals(end)])
xlim([D.T(1)-.01 D.T(end)])

axis xy


