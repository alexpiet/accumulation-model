# Accumulation model
Hello!

This code package implements an analytical solution to a drift-diffusion-advection
(linear feedback) model for pulse-based evidence. It can be used to do model
fitting via maximum likelihood estimation. In addition, it can compute the
so-called backwards pass distribution which gives a better estimation of the
latent variable at each time point by using the final choice of the agent.

For questions, please contact alexpiet@gmail.com. This code is being provided
for the community, and is not documented as well as it could be. Sorry!

### Getting started with sample data
I've attached one sample datafile "H033.mat" which is a rat from Piet et al,
Nature Communications, 2018. If you want to fit an agent's behavior, start with
the file "fit_rat_analytical.m". If you want to use the backwards pass, start with
"dev_script.m" and "accumulation_model.m". The file "backwards_pass_solution.pdf"
documents how the code works.

### Citation
If this code is useful to you, you can cite Piet, AT, El Hady A, & Brody CD. Rats adopt the optimal
timescale for evidence integration in a dynamic environment. Nature Communications, 2018.
That paper used the forward model to model rat behavior in a decision making task.
A manuscript in progress uses the backwards pass.
