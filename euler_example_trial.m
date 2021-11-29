% set agent parameters
param_names = {'\lambda', '\sigma_a', '\sigma_s', '\sigma_i', '\phi',  '\tau', 'bias', 'lapse'};
params = [0 1 1 0 .5 .1 0 0];
lambda = params(1); var_a = params(2); var_s = params(3); var_init = params(4);
phi = params(5); tau_phi = params(6); bias = params(7); lapse = params(8); 
bound = Inf; 
n_particles = 100;
n = n_particles;
n_to_plot = 5;

% make clicks
total_rate = 40;
gamma = .5;
l_rate = total_rate ./ (exp( gamma) + 1);
r_rate = total_rate - l_rate;

T  = 1;
dt = 1e-3;
tvec = 0:dt:T;
% set up the random seed that determines the noise realizations
rng(2)
lbupvec = rand(size(tvec)) < l_rate*dt;
rbupvec = rand(size(tvec)) < r_rate*dt;
lb = tvec(find(lbupvec));
rb = tvec(find(rbupvec));

% apply click adaptation
[cl, cr]    = make_adapted_cat_clicks(lb, rb, phi, tau_phi);
[difflr, sumlr] = make_click_inputs35(tvec, lb, rb, cl, cr);

% simulate particles
init_noise = sqrt(var_init).*randn(n,1);
a = zeros(n,length(tvec));
a(:,1) = init_noise;
for ii=2:length(tvec)
    last_a = a(:,ii-1);
    to_update = abs(last_a ) < bound;
    n = sum(to_update);
    new_a = last_a(to_update) ...
        + dt*lambda.*last_a(to_update) ...
        + difflr(ii-1) + sqrt(sumlr(ii-1)*var_s).*randn(n,1) ...
        + sqrt(var_a*dt).*randn(n,1);
    hit_bound = abs(new_a ) >= bound;
    new_a(hit_bound) = bound*sign(new_a(hit_bound));
    a(to_update,ii) = new_a;
    a(~to_update,ii) = a(~to_update,ii-1);
            
end

% plot the clicks and the particles
fh = figure(1); clf
set(fh,'units','inches','position',[5 5 10 5])
subplot(3,3,[1 2])
left_color = [.75 1 .75];
right_color = [.75 .75 1];
plot(repmat(rb,2,1),[0; 1], '-','color',left_color,'linewidth',1.5)
hold on
plot(repmat(lb,2,1),[0; -1], '-','color',right_color,'linewidth',1.5)
text(-.05, 1.5, 'clicks','fontsize',15)
text(.1, 1.5, 'left','fontsize',15, 'color', left_color)
text(.4, 1.5, 'right','fontsize',15, 'color', right_color)

left_color = [.05 .85 .05];
right_color = [.05 .05 .85];
plot(repmat(rb,2,1),[zeros(size(cr)); cr], '-','color',left_color,'linewidth',2)
hold on
plot(repmat(lb,2,1),[zeros(size(cl)); -cl], '-','color',right_color,'linewidth',2)
axis off
text(.2, 1.5, 'adapted','fontsize',15, 'color', left_color)
text(.5, 1.5, 'adapted','fontsize',15, 'color', right_color)
box off

ax2 = subplot(3,3,[4 5 7 8])
% plot particles
plot(tvec, a(1:n_to_plot,:) ,'linewidth', 1.5)
hold on
plot(tvec([1 end]), bias*[1 1], '--k')
box off
ylabel('a')
xlabel('time (s)')
final_a = a(:,end);
a_greater_than_bias = final_a > bias;

ylim(1.25*ylim)
text(0.05, .8*max(ylim), sprintf(['\\lambda=%.1f, \\sigma^2_a=%.1f, \\sigma^2_i=%.1f, '...
    '\\sigma^2_s=%.1f, \\phi=%.1f, \\tau_{\\phi}=%.2f, bias=%.1f'],...
    lambda, var_a, var_init, var_s, phi, tau_phi, bias),'fontsize',13)

text(.1, 1.0*max(ylim), sprintf('a>bias for %i/%i particles with parameters', sum(a_greater_than_bias), n_particles), 'fontsize', 15)


buptimes = [lb rb];
streamIdx = [-ones(size(lb)) ones(size(rb))]
[buptimes,nantimes,streamIdx] = vectorize_clicks({lb},{rb})
[~, ma, va, ~, ~, pr] = compute_LL_vectorized(buptimes,streamIdx,...
            T, 1, params(1:8), 'nantimes', nantimes)
        
subplot(3,3,[6 9])
h = histogram(final_a,10,'normalization','pdf','facecolor',[1 1 1].*.75)

hold on
xx = linspace(min(final_a),max(final_a),100);
plot(xx,normpdf(xx,ma,sqrt(va)),'r-')
view(-90,90);
box off
xlim(get(ax2,'ylim'))
title('final value of a','fontweight','normal')
%xlabel('a')
ylabel('% realizations')
axis ij
plot([1 1]*bias,ylim,'k--')
text(bias,max(ylim)*.8, 'bias','fontsize',12)
legend('simulation','model','location','northeast')
