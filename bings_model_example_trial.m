%% user parameters to modify
% click parameters
total_rate = 40; % left and right rates sum to this value
gamma = .5; % log ratio of the left and right gamma = log(r_R/r_L)
trial_duration  = 1;
dt = 1e-3;

% agent parameters
lambda = -.4; % positive values -> instability; negative values -> leak
var_a = 0; % variance of noise applied at each time step
var_s = 36.6/total_rate; % variance of noise applied with each click input - multiply by total_rate to get bing's units
var_init = .1; % variance of noise applied at beginning of trial
phi = .15; % adaptation strength: >1 -> facilitation; <1-> depression
tau_phi = .14; % timescale of adaptation
bias = -.18; % go right if a > bias
lapse = .1; % fraction of time where agent chooses randomly; not used here.
bound = 5; % height of absorbing bound; will not be reflected by analytical model

w_l = 1;
w_r = 1;

% simulation parameters
n_particles = 500; % number of particles to simulate
n_to_plot = 5; % number of particles to plot

% don't modify beyond here 
params = [lambda, var_a, var_s, var_init, phi, tau_phi, bias, lapse];

% make clicks
l_rate = total_rate ./ (exp( gamma) + 1);
r_rate = total_rate - l_rate;
tvec = 0:dt:trial_duration;

% set up the random seed that determines the noise realizations
rng(2)
lbupvec = rand(size(tvec)) < l_rate*dt;
rbupvec = rand(size(tvec)) < r_rate*dt;
lbupvec(1) = 1;
rbupvec(1) = 1;
lb = tvec(find(lbupvec));
rb = tvec(find(rbupvec));

% apply click adaptation
[cl, cr]    = make_adapted_cat_clicks(lb, rb, phi, tau_phi);
wcl = cl*w_l;
wcr = cr*w_r;
[difflr, sumlr] =  make_click_inputs35(tvec, lb, rb, wcl, wcr);

% simulate particles
n = n_particles;
init_noise = sqrt(var_init).*randn(n,1);
a = zeros(n,length(tvec));
a(:,1) = init_noise;
for ii=2:length(tvec)
    last_a = a(:,ii-1);
    to_update = abs(last_a ) < bound;
    n = sum(to_update);
    new_a = last_a(to_update) ... % previous value
        + dt*lambda.*last_a(to_update) ... % effect of leak/instability
        + difflr(ii-1) ... % effect of clicks in this timestep
        + sqrt(sumlr(ii-1)*var_s).*randn(n,1) ... % noise due to clicks
        + sqrt(var_a*dt).*randn(n,1); % accumulator noise 
    hit_bound = abs(new_a ) >= bound;
    new_a(hit_bound) = bound*sign(new_a(hit_bound));
    a(to_update,ii) = new_a;
    a(~to_update,ii) = a(~to_update,ii-1);
            
end

% plot the clicks and the particles
fh = figure(1); clf
set(fh,'units','inches','position',[0 5 12 6])
subplot(3,3,[1 2])
right_color = [.05 .85 .05];
left_color = [.05 .05 .85];
lighter = @(x) hsv2rgb(rgb2hsv(x) .* [1 .2 1]);
right_color_light = lighter(right_color);
left_color_light = lighter(left_color);
if ~isempty(rb)
    plot(repmat(rb,2,1),[0; 1], '-','color',...
    right_color_light,'linewidth',1.5)
end
hold on
if ~isempty(lb)
    plot(repmat(lb,2,1),[0; -1], '-','color', left_color_light,'linewidth',1.5)
end
text(-.05, 1.5, 'clicks','fontsize',15)
text(.1, 1.5, 'left','fontsize',15, 'color', left_color_light)
text(.4, 1.5, 'right','fontsize',15, 'color', right_color_light)


plot(repmat(rb,2,1),[zeros(size(cr)); cr], '-','color', right_color,'linewidth',2)
hold on
plot(repmat(lb,2,1),[zeros(size(cl)); -cl], '-','color', left_color,'linewidth',2)
axis off
text(.2, 1.5, 'adapted','fontsize',15, 'color', left_color)
text(.5, 1.5, 'adapted','fontsize',15, 'color', right_color)
box off

text(1.125, 1.2, ['$\log \frac{r_R}{r_L}=$' num2str(gamma)] ,...
    'fontsize',17, 'interpreter','latex')
text(1.125, .5, sprintf('$r_R=%.1f$ clicks $s^{-1}$',r_rate), ...
    'fontsize',17, 'interpreter','latex','color',right_color)
text(1.125, -.5, sprintf('$r_L=%.1f$ clicks $s^{-1}$',l_rate), ...
    'fontsize',17, 'interpreter','latex','color',left_color)


% plot particles

ax2 = subplot(3,3,[4 5 7 8])

plot(tvec, a(1:n_to_plot,:) ,'linewidth', 1.5,'color',[1 1 1].*.4)
hold on
plot(tvec([1 end]), bias*[1 1], '--k')
box off
ylabel('accumulation value, a','fontsize',18)
xlabel('time (s)')
final_a = a(:,end);
a_greater_than_bias = final_a > bias;

ylim([min(final_a) max(final_a)])
ylims = ylim(ax2)
if ylims(2) > .9*bound
    ylims(2) = 1.3*bound
end
ylim(ax2,ylims)

text(0.05, .9*max(ylim), sprintf(['\\lambda=%.1f, \\sigma^2_a=%.1f, \\sigma^2_i=%.1f, '...
    '\\sigma^2_s=%.1f, \\phi=%.1f, \\tau_{\\phi}=%.2f, bias=%.1f'],...
    lambda, var_a, var_init, var_s, phi, tau_phi, bias),'fontsize',13)

text(.1, 1.1*max(ylim), sprintf('a>bias for %i/%i particles with parameters', sum(a_greater_than_bias), n_particles), 'fontsize', 15)


buptimes = [lb rb];
streamIdx = [-ones(size(lb)) ones(size(rb))]
[buptimes,nantimes,streamIdx] = vectorize_clicks({lb},{rb})
[~, ma, va, ~, ~, pr] = compute_LL_vectorized(buptimes,streamIdx,...
            trial_duration, 1, params(1:8), 'nantimes', nantimes)
        
subplot(3,3,[6 9])
h = histogram(final_a,10,'normalization','pdf','facecolor',[1 1 1].*.75)

hold on
xx = linspace(min(final_a),max(final_a),100);
plot(xx,normpdf(xx,ma,sqrt(va)),'r-')
view(-90,90);
box off
xlim(get(ax2,'ylim'))
title('Final value of a','fontweight','normal')
%xlabel('a')
ylabel('% realizations')
axis ij
plot([1 1]*bias,ylim,'k--')
text(bias,max(ylim)*.8, 'bias','fontsize',12)
hl = legend('simulation (has bound)','model (no bound)','location','southeast')
box(hl,'off')
hl.Position = hl.Position + [.05 0 0 0]


