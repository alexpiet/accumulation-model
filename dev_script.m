%close all; clear all;

% load example fit and dataset
dp = set_dyn_path(1);
ex_rat = 'H084';

fit = fit_rat_analytical(ex_rat,'data_dir',dp.data_dir,'results_dir',dp.model_fits_dir);
load(fullfile(dp.data_dir,ex_rat));
%%

% set up parameters
TN = 1;
data(TN).pokedR  = 0; % Which choice to compute posterior for
params = fit.final;
params(1) = 0;
%params(2) = 1e-4; % accumulation noise
%params(3) = 1e-4; % sensory noise
%params(4) = 1e-4; % initial noise
params(7) = 0;
p.compute_full  = 1; 
p.compute_back  = 1;
p.compute_delta = 5;
p.save_new_ref  = 0;    
p.dt            = 1e-4;
%p.check_final   = 1;
p.b             = 50;
p.da            = 1;
p.avals         = -p.b:p.da:p.b;
p.a_sign = 2 * data(TN).pokedR - 1; % added by TB
p.da_grid       = 1; % uncommented by TB
p.a_grid        = -p.b:p.da_grid:p.b + params(7); % uncommented by TB
%%p.a_grid        = [-p.b:p.da_grid:-2 -1.75:0.25:-1 0:p.da_grid:p.b] + params(7);
%%p.a_grid        = [-p.b:2:-10 -9:p.da_grid:p.b] + params(7);
%p.a_grid_s = [p.da_grid (diff(p.a_grid(1:end-1)) + diff(p.a_grid(2:end)))/2 p.da_grid];
%% da grid notes
%% 0 bin handled elsewhere
%% fix grid sizes in compute_particles (3 places)
%% set scale factors in compute_backwards_pass
p.d_dex         = 6 * p.a_sign; % abs value of delta-mode to compute
p.n             = 10000;
p.error_tolerance = 0.0001;
computet        = 0;
plott           = 0;
particlet       = 0;

if length(unique(diff(p.avals))) > 1
    error('p.avals is not a uniform grid! I cannot verify solution for a non uniform evaluation grid. Proceed with EXTREME caution!');
end
if p.da > 1
    error('p.da > 1 may lead to unstable behavior at some point. Proceed with caution!');
end

%if p.da_grid > p.da
%    error('delta mode spacing larger than evaluation spacing, will create funky results. Proceed with caution!');
%end

if params(3) < 0.5
    error('WARNING, you might be hitting the low-variance case. See documentation for details, and for advice on how to proceed.')
end
% compute foward pass distribution
tic
[cl, cr]    = make_adapted_cat_clicks(data(TN).leftbups, data(TN).rightbups, params(5), params(6));
data(TN).clicks      = [-cl  +cr];
data(TN).times       = [data(TN).leftbups data(TN).rightbups];
data(TN).times       = round(data(TN).times/p.dt)*p.dt;
data(TN).T           = round(data(TN).T/p.dt)*p.dt;
data(TN).numsteps    = round(data(TN).T/p.dt);
data(TN).dtimes      = round(data(TN).times*(1/p.dt));
[forward] = compute_full_trial(data(TN),params,p); 
computet  = computet +toc;

% Visualize forward pass
tic
forward = compute_pdf(forward,p.avals,p);
forward.plot_zero_line = true;
forward.plot_mean_line = false;
forward.title = 'Forward';
forward.right_click_marker = '|';
forward.left_click_marker = '|';
forward.left_clicks = data(TN).leftbups;
forward.right_clicks = data(TN).rightbups;
forward.right_click_y = 25;
forward.left_click_y = -25;
forward.click_height = 5;
%%
[~, alimi] = max(abs(p.a_grid));
alim = p.a_grid(alimi);
fh = figure(1); 
clf; 
subplot(2,4,1);
plot_pdf(forward);

plott = plott + toc;
% compute backwards pass
if p.compute_back;
tic
p = get_grid(data(TN), forward, params,p);
[back, posterior] = compute_backwards_pass(data(TN),params,p,forward);
computet = computet + toc;

%%
bdex = find(back.a_grid == p.d_dex); % which mode to compute
back = compute_pdf(back,p.avals,p,'mode',bdex);
posterior = compute_pdf(posterior, p.avals,p,'mode',bdex);
back.s = ones(size(back.ma))./50; %% DOES NOTHING
backF = compute_pdf(back,p.avals,p,'mixture');
posteriorF = compute_pdf(posterior, p.avals,p,'mixture');
particle = compute_particles(data(TN), params, p);
%%
% visualize backwards pass - one delta solution
tic
back.plot_zero_line = true;
back.plot_mean_line = false;
back.title = sprintf('Backward \\delta, a_N = %i',p.d_dex);
back.right_click_marker = '|';
back.left_click_marker = '|';
back.left_clicks = data(TN).leftbups;
back.right_clicks = data(TN).rightbups;
back.right_click_y = 25;
back.left_click_y = -25;
back.click_height = 5;
subplot(2,4,2);
plot_pdf(back);

% visualize posterior - one delta solution
subplot(2,4,3);
posterior.plot_zero_line = true;
posterior.plot_mean_line = false;
posterior.title = sprintf('Posterior \\delta, a_N = %i',p.d_dex);
posterior.right_click_marker = '|';
posterior.left_click_marker = '|';
posterior.left_clicks = data(TN).leftbups;
posterior.right_clicks = data(TN).rightbups;
posterior.right_click_y = 25;
posterior.left_click_y = -25;
posterior.click_height = 5;
plot_pdf(posterior);
% 
% % visualize entire backwards pass
% backF.plot_zero_line = true;
% backF.plot_mean_line = false;
% backF.title = 'Backward Full';
% backF.right_click_marker = '|';
% backF.left_click_marker = '|';
% backF.left_clicks = data(TN).leftbups;
% backF.right_clicks = data(TN).rightbups;
% backF.right_click_y = 25;
% backF.left_click_y = -25;
% backF.click_height = 5;
% subplot(2,4,5);
% plot_pdf(backF);

% visualize entire posterior 
subplot(2,4,4);
plott = plott + toc;
tic
computet = computet + toc;

tic
posteriorF.plot_zero_line = true;
posteriorF.plot_mean_line = false;
if p.a_sign
    posteriorF.title = sprintf('Posterior Full, a_N = U(0,%i)',alim);
else
    posteriorF.title = sprintf('Posterior Full, a_N = U(%i,0)',alim);
end
posteriorF.right_click_marker = '|';
posteriorF.left_click_marker = '|';
posteriorF.left_clicks = data(TN).leftbups;
posteriorF.right_clicks = data(TN).rightbups;
posteriorF.right_click_y = 25;
posteriorF.left_click_y = -25;
posteriorF.click_height = 5;
plot_pdf(posteriorF);
plott = plott + toc;
end

%% particles
tic
particlet = particlet + toc;
tic
subplot(2,4,5);
forwardP.pdf = particle.fpdf;
forwardP.avals = particle.avals;
forwardP.T = particle.T;
forwardP.plot_zero_line = true;
forwardP.plot_mean_line = false;
forwardP.title = 'Particle Forward';
forwardP.right_click_marker = '|';
forwardP.left_click_marker = '|';
forwardP.left_clicks = data(TN).leftbups;
forwardP.right_clicks = data(TN).rightbups;
forwardP.right_click_y = 25;
forwardP.left_click_y = -25;
forwardP.click_height = 5;
plot_pdf(forwardP)

subplot(2,4,6);
backP.pdf = particle.bpdf;
backP.avals = particle.avals;
backP.T = particle.T;
backP.plot_zero_line = true;
backP.plot_mean_line = false;
%backP.title = {'Particle Backward \delta'};
backP.title = sprintf('Particle Backward \\delta, a_N = %i',p.d_dex);
backP.right_click_marker = '|';
backP.left_click_marker = '|';
backP.left_clicks = data(TN).leftbups;
backP.right_clicks = data(TN).rightbups;
backP.right_click_y = 25;
backP.left_click_y = -25;
backP.click_height = 5;
plot_pdf(backP)

subplot(2,4,7);
posteriorD.pdf = particle.dpdf;
posteriorD.avals = particle.avals;
posteriorD.T = particle.T;
posteriorD.plot_zero_line = true;
posteriorD.plot_mean_line = false;
posteriorD.title = sprintf('Particle Posterior \\delta, a_N = %i',p.d_dex);
posteriorD.right_click_marker = '|';
posteriorD.left_click_marker = '|';
posteriorD.left_clicks = data(TN).leftbups;
posteriorD.right_clicks = data(TN).rightbups;
posteriorD.right_click_y = 25;
posteriorD.left_click_y = -25;
posteriorD.click_height = 5;
plot_pdf(posteriorD)


ax = subplot(2,4,8);
posteriorP.pdf = particle.ppdf;
posteriorP.avals = particle.avals;
posteriorP.T = particle.T;
posteriorP.plot_zero_line = true;
posteriorP.plot_mean_line = false;

if p.a_sign
    posteriorP.title = sprintf('Posterior Particle Full, a_N = U(0,%i)',alim);
else
    posteriorP.title = sprintf('Posterior Particle Full, a_N = U(%i,0)',alim);
end
posteriorP.right_click_marker = '|';
posteriorP.left_click_marker = '|';
posteriorP.left_clicks = data(TN).leftbups;
posteriorP.right_clicks = data(TN).rightbups;
posteriorP.right_click_y = 25;
posteriorP.left_click_y = -25;
posteriorP.click_height = 5;
plot_pdf(posteriorP)
plott = plott + toc;
axpos = get(ax,'position');
colorbar
drawnow
set(ax,'position',axpos)
%%
allax = findall(fh,'type','axes');
set(allax,'Box','off');
ylim(allax,[-40 40]);
end_color = dp.model_color;
colormap(colormapLinear(end_color).^2)
for aa = 1:length(allax)
    caxis(allax(aa),[0 .1])
    view(allax(aa),[90 90])
end

fht = 2.5 * 2.5;
fw  = 4/3*2.5 * 3.75;
set(fh, 'position', [0 10 fw fht],'papersize',[fw fht],'paperpositionmode','auto')
if p.a_sign > 0
    print(fullfile(dp.fig_dir, 'dist_check_supp_right'), '-dsvg','-painters')
else
    print(fullfile(dp.fig_dir, 'dist_check_supp_left'), '-dsvg','-painters')
end
%%
keyboard

%% Numerical checks FORWARD
figure;
pma = mean(particle.a);
diffma = sum(abs(pma - forward.ma))./length(forward.T);
subplot(3,3,1);hold on;
plot(forward.T, forward.ma,'b');
plot(forward.T, pma,'r')
title('forward mean'); xlabel('Time (s)')

pause(.1)
pva = var(particle.a);
diffva = sum(abs(pva - forward.va))./length(forward.T);
subplot(3,3,2);hold on;
plot(forward.T, forward.va,'b');
plot(forward.T, pva,'r')
title('forward variance'); xlabel('Time (s)')

subplot(3,3,3); hold on;
plot(p.avals, forward.pdf(1000,:), 'b')
plot(p.avals, forwardP.pdf(1000,:), 'r')
title('forward slice at 0.1s'); xlabel(' a value')

subplot(3,3,4); hold on;
plot(p.avals, forward.pdf(3000,:), 'b')
plot(p.avals, forwardP.pdf(3000,:), 'r')
title('forward slice at 0.3s'); xlabel(' a value')

subplot(3,3,5); hold on;
plot(p.avals, forward.pdf(5000,:), 'b')
plot(p.avals, forwardP.pdf(5000,:), 'r')
title('forward slice at 0.5s'); xlabel(' a value')

subplot(3,3,6); hold on;
plot(p.avals, forward.pdf(7000,:), 'b')
plot(p.avals, forwardP.pdf(7000,:), 'r')
title('forward slice at 0.7s'); xlabel(' a value')

subplot(3,3,7); hold on;
plot(p.avals, forward.pdf(9000,:), 'b')
plot(p.avals, forwardP.pdf(9000,:), 'r')
title('forward slice at 0.9s'); xlabel(' a value')

subplot(3,3,8); hold on;
plot(p.avals, forward.pdf(end-1,:), 'b')
plot(p.avals, forwardP.pdf(end-1,:), 'r')
title('forward slice at end-1'); xlabel(' a value')

subplot(3,3,9); hold on;
plot(p.avals, forward.pdf(end,:), 'b')
plot(p.avals, forwardP.pdf(end,:), 'r')
title('forward slice at end'); xlabel(' a value')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 9];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [9 9 ];
%print([pwd '/dev/figs/forward_check'], '-dpdf')

%% Numerical checks BACKWARDS-DELTA
figure;
bdma = mean(particle.b(:,:));
subplot(3,3,1);hold on;
plot(forward.T, back.ma(bdex,:),'b');
plot(forward.T, bdma,'r')
title('backward-d mean'); xlabel('Time (s)')

pause(.1)
bdva = var(particle.b(:,:));
subplot(3,3,2);hold on;
plot(forward.T, back.va(bdex,:),'b');
plot(forward.T, bdva,'r')
title('backward-d variance'); xlabel('Time (s)')

subplot(3,3,3); hold on;
plot(p.avals, back.pdf(1000,:), 'b')
plot(p.avals, backP.pdf(1000,:), 'r')
title('backward-d slice at 0.1s'); xlabel(' a value')

subplot(3,3,4); hold on;
plot(p.avals, back.pdf(3000,:), 'b')
plot(p.avals, backP.pdf(3000,:), 'r')
title('backward-d slice at 0.3s'); xlabel(' a value')

subplot(3,3,5); hold on;
plot(p.avals, back.pdf(5000,:), 'b')
plot(p.avals, backP.pdf(5000,:), 'r')
title('backward-d slice at 0.5s'); xlabel(' a value')

subplot(3,3,6); hold on;
plot(p.avals, back.pdf(7000,:), 'b')
plot(p.avals, backP.pdf(7000,:), 'r')
title('backward-d slice at 0.7s'); xlabel(' a value')

subplot(3,3,7); hold on;
plot(p.avals, back.pdf(8000,:), 'b')
plot(p.avals, backP.pdf(8000,:), 'r')
title('backward-d slice at 0.8s'); xlabel(' a value')

subplot(3,3,8); hold on;
plot(p.avals, back.pdf(9000,:), 'b')
plot(p.avals, backP.pdf(9000,:), 'r')
title('backward-d slice at 0.9s'); xlabel(' a value')

subplot(3,3,9); hold on;
plot(p.avals, back.pdf(end-1,:), 'b')
plot(p.avals, backP.pdf(end-1,:), 'r')
title('backward-d slice at end'); xlabel(' a value')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 9];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [9 9 ];
%print([pwd '/dev/figs/backward_delta_check'], '-dpdf')

%% Numerical checks  POSTERIOR-DELTA
figure;
pdma = mean(particle.a(particle.dDex,:));
diffdma = sum(abs(pdma - posterior.ma(bdex,:)))./length(forward.T);
subplot(3,3,1);hold on;
plot(forward.T, posterior.ma(bdex,:),'b');
plot(forward.T, pdma,'r')
title('posterior-d mean'); xlabel('Time (s)')

pause(.1)
pdva = var(particle.a(particle.dDex,:));
diffdva = sum(abs(pdva - posterior.va(bdex,:)))./length(forward.T);
subplot(3,3,2);hold on;
plot(forward.T, posterior.va(bdex,:),'b');
plot(forward.T, pdva,'r')
title('posterior-d variance'); xlabel('Time (s)')

subplot(3,3,3); hold on;
plot(p.avals, posterior.pdf(1000,:), 'b')
plot(p.avals, posteriorD.pdf(1000,:), 'r')
title('posterior-d slice at 0.1s'); xlabel(' a value')

subplot(3,3,4); hold on;
plot(p.avals, posterior.pdf(3000,:), 'b')
plot(p.avals, posteriorD.pdf(3000,:), 'r')
title('posterior-d slice at 0.3s'); xlabel(' a value')

subplot(3,3,5); hold on;
plot(p.avals, posterior.pdf(5000,:), 'b')
plot(p.avals, posteriorD.pdf(5000,:), 'r')
title('posterior-d slice at 0.5s'); xlabel(' a value')

subplot(3,3,6); hold on;
plot(p.avals, posterior.pdf(7000,:), 'b')
plot(p.avals, posteriorD.pdf(7000,:), 'r')
title('posterior-d slice at 0.7s'); xlabel(' a value')

subplot(3,3,7); hold on;
plot(p.avals, posterior.pdf(8000,:), 'b')
plot(p.avals, posteriorD.pdf(8000,:), 'r')
title('posterior-d slice at 0.8s'); xlabel(' a value')

subplot(3,3,8); hold on;
plot(p.avals, posterior.pdf(9000,:), 'b')
plot(p.avals, posteriorD.pdf(9000,:), 'r')
title('posterior-d slice at 0.9s'); xlabel(' a value')

subplot(3,3,9); hold on;
plot(p.avals, posterior.pdf(end-1,:), 'b')
plot(p.avals, posteriorD.pdf(end-1,:), 'r')
title('posterior-d slice at end'); xlabel(' a value')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 9];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [9 9 ];
%print([pwd '/dev/figs/posterior_delta_check'], '-dpdf')

% POSTERIOR-FULL
% load reference posterior
if TN == 1
load(fullfile(dp.data_dir,'reference_posterior.mat'));
end

figure;
subplot(3,3,1); hold on;
for i=1:length(posteriorF.M)
    maF(i) = sum(p.avals'.*posteriorF.PDF{i}.*p.da);
end
posteriorF.maF = maF;
plot(forward.T, posteriorF.maF,'b')
pfma = mean(particle.a(particle.DEX,:));
plot(forward.T, pfma, 'r')
if TN == 1
plot(forward.T, ref.pfma, 'r--')
plot(forward.T, ref.posteriorF.maF,'b--')
end
xlabel('time (s)')
title('posterior mean');

subplot(3,3,2); hold on;
for i=1:length(posteriorF.M)
    fx = (p.avals' - posteriorF.maF(i)).^2;
    vaF(i) = sum(fx.*posteriorF.PDF{i}.*p.da);
end
posteriorF.vaF = vaF;
plot(forward.T, posteriorF.vaF,'b')
pfva = var(particle.a(particle.DEX,:));
plot(forward.T, pfva, 'r')
if TN == 1
plot(forward.T, ref.pfva, 'r--')
plot(forward.T, ref.posteriorF.vaF,'b--')
end
xlabel('time (s)')
title('posterior variance')

subplot(3,3,3); hold on;
plot(p.avals, posteriorF.pdf(1000,:), 'b')
plot(p.avals, posteriorP.pdf(1000,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(1000,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(1000,:), 'r--')
end

xlabel(' a value')
title('posterior slice at 0.1s')

subplot(3,3,4); hold on;
plot(p.avals, posteriorF.pdf(3000,:), 'b')
plot(p.avals, posteriorP.pdf(3000,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(3000,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(3000,:), 'r--')
end
xlabel(' a value')
title('posterior slice at 0.3s')

subplot(3,3,5); hold on;
plot(p.avals, posteriorF.pdf(5000,:), 'b')
plot(p.avals, posteriorP.pdf(5000,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(5000,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(5000,:), 'r--')
end
xlabel(' a value')
title('posterior slice at 0.5s')

subplot(3,3,6); hold on;
plot(p.avals, posteriorF.pdf(7000,:), 'b')
plot(p.avals, posteriorP.pdf(7000,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(7000,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(7000,:), 'r--')
end
xlabel(' a value')
title('posterior slice at 0.7s')

subplot(3,3,7); hold on;
plot(p.avals, posteriorF.pdf(9000,:), 'b')
plot(p.avals, posteriorP.pdf(9000,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(9000,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(9000,:), 'r--')
end
xlabel(' a value')
title('posterior slice at 0.9s')

subplot(3,3,8); hold on;
plot(p.avals, posteriorF.pdf(end-1,:), 'b')
plot(p.avals, posteriorP.pdf(end-1,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(end-1,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(end-1,:), 'r--')
end
xlabel(' a value')
title('posterior slice at end-1')

subplot(3,3,9); hold on;
plot(p.avals, posteriorF.pdf(end,:), 'b')
plot(p.avals, posteriorP.pdf(end,:), 'r')
if TN == 1
plot(ref.posteriorF.avals, ref.posteriorF.pdf(end,:), 'b--')
plot(ref.posteriorF.avals, ref.posteriorP.pdf(end,:), 'r--')
end
xlabel(' a value')
title('posterior slice at end')

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 9];
fig.PaperPositionMode = 'Manual';
fig.PaperSize = [9 9 ];
%print([pwd '/dev/figs/posterior_check'], '-dpdf')

disp(['Time:'])
disp(['Analytic Computation : ' num2str(computet) ' seconds'])
disp(['Particle Simulations : ' num2str(particlet) ' seconds'])
disp(['Dist. Visualization  : ' num2str(plott) ' seconds'])


if p.save_new_ref;
    ref.pfma = pfma;
    ref.pfva = pfva;
    ref.posteriorF = posteriorF;
    ref.posteriorP = posteriorP;
    save(fullfile(dp.data_dir,'reference_posterior.mat'),'ref');
end   

