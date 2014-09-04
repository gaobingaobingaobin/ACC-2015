clear all 
close all 
clc

% verify that required toolboxes are installed 
check_system_requirements(); 

% set colors for plots 
berkeley_colors = ...
 1/256*[ 45,  99, 127; 
        224, 158,  25; 
          0,   0,   0;
        194, 185, 167;
        217, 102, 31;
        185, 211, 182]; 

    
%% Specify system model 

% initialize model object 
model = linear_exchange_model; 

% define model parameters
syms R1P R1L R1A kPL kPA kTRANS 
% define input parameters 
syms t0 alpha_1 beta_1 A0 
% define noise parameters 
syms sigma_1 sigma_2 sigma_3
% define initial state parameters
syms P0 L0 

% parameters of interest 
% (those for which we wish to compute an estimate with minimal variance) 
model.parameters_of_interest = [kPL kTRANS]; 
model.parameters_of_interest_nominal_values = [0.02 0.04]; 

% nuisance parameters
% (those parameters that are unknown but whose estimates we only care about
% insofar as they allow us to estamate the parameters of interest) 
model.nuisance_parameters = [alpha_1 beta_1 A0];
model.nuisance_parameters_nominal_values = [ 2  5  1]; 

% known parameters
% (those whose values are assumed to be known constants) 
model.known_parameters = [R1P R1L R1A t0 P0 L0]; 
model.known_parameter_values = [1/35 1/30 1/25 0 0 0];  

% define system matrices for differential eq. dx/dt = A*x(t) + B*u(t)
% the input should be the same chemical compound as the first state 

% two-site exchange model 
model.A = [ -kPL-R1P  0   ;
             kPL     -R1L];   
model.B = [kTRANS; 0]; 

% three-site exhange model 
% model.A = [ -kPL-kPA-R1P  0    0;
%             kPL         -R1L  0;
%             kPA          0   -R1A];   
% model.B = [kTRANS; 0; 0]; 

% define input function shape  
model.u = @(t) A0 * (t - t0)^alpha_1 *exp(-(t - t0)/beta_1); % gamma-variate input  
% model.u = @(t) 10*rectangularPulse(0, 15, t);              % boxcar input 

% define initial condition 
model.x0 = [P0; L0]; 
% model.x0 = [P0; L0; 0]; 

% define repetition time
model.TR = 2; 

% define number of acquisitions 
model.N = 25; 

% choose noise type
model.noise_type = 'Rician';
% model.noise_type = 'None';

% choose noise magnitude  
model.noise_parameters = [0.01 0.01 0.01]; 
% model.noise_parameters = [0.01 0.01 0.01 0.1]; 

% choose flip angle input matrix 
%   This allows you to set linear equality constraints on the flip angles
%   for example setting: 
%
%      model.flip_angle_input_matrix = [1 0; 
%                                       0 1; 
%                                       1 0]; 
%
%   fixes the first and third flip angles to be equal one another. 
%   The attribute model.flip_angle_input_matrix must have m+n rows. 
%   Consider defaulting to
% 
%      model.flip_angle_input_matrix = eye(model.m + model.n) 
% 
%   if you wish to compute all flip angles separately. 

model.flip_angle_input_matrix = eye(model.m + model.n);                               

% model.flip_angle_input_matrix = [1 0; 
%                                 0 1; 
%                                 1 0]; 
% model.flip_angle_input_matrix = eye(model.m + model.n)                              

% choose design criterion 
design_criterion = 'D-optimal'; 
% design_criterion = 'E-optimal'; 
% design_criterion = 'A-optimal'; 
% design_criterion = 'T-optimal'; 
% design_criterion = 'totalSNR'; 

% discretize model (doing this in advance makes things run faster) 
model = discretize(model);  

% compute sensitivities (doing this in advance makes things run faster)
if ~model.sensitivities_computed ...
        && (strcmp(design_criterion, 'D-optimal') ...
            || strcmp(design_criterion, 'E-optimal') ...
            || strcmp(design_criterion, 'A-optimal') ...
            || strcmp(design_criterion, 'T-optimal') ...
           )
    model = sensitivities(model);  
end

    
%% Figure 1: Numerically computed integral phi 

% compute the function phi 
phi = compute_phi(); 
z_vals = linspace(0, 4, 500); 

% plot phi 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(z_vals, phi(z_vals), 'LineWidth', 2); 
axis([0 4 0 1]); 
xlabel('x/\sigma')
ylabel('\phi')
title('Numerically computed integral \phi for various SNR values') 
print(h, '-dpdf', 'fig1.pdf')


%% Figure 2: Constant flip angle design for D- E- and A-optimal objectives

% define objective functions 
f_D = @(theta) log(abs(det(fisher_information(theta*ones(model.N, model.n + model.m), model, phi))));
f_E = @(theta) max(eig(fisher_information(theta*ones(model.N, model.n + model.m), model, phi)));
f_A = @(theta) 1/trace(inv(fisher_information(theta*ones(model.N, model.n + model.m), model, phi)));

% plot objective function values 
theta_vals = linspace(0.01, pi/2, 181); 

% compute objective function values
fD_vals = zeros(size(theta_vals)); 
fE_vals = zeros(size(theta_vals)); 
fA_vals = zeros(size(theta_vals)); 
for i=1:length(fD_vals)
    fD_vals(i) = f_D(theta_vals(i)); 
    fE_vals(i) = f_E(theta_vals(i)); 
    fA_vals(i) = f_A(theta_vals(i)); 
end

% plot D-optimal objective 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(theta_vals, fD_vals, 'LineWidth', 2); 
axis([0 pi/2 0 25]); 
xlabel('\theta')
ylabel('log(det(I))')
title('D-optimal objective function')
print(h, '-dpdf', 'fig2a.pdf')

% plot E-optimal objective 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(theta_vals, fE_vals, 'LineWidth', 2); 
axis([0 pi/2 0 8e05]); 
xlabel('\theta')
ylabel('max(eig(I))')
title('E-optimal objective function')
print(h, '-dpdf', 'fig2b.pdf')

% plot A-optimal objective 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(theta_vals, fA_vals, 'LineWidth', 2); 
axis([0 pi/2 0 3e04]); 
xlabel('\theta')
ylabel('1/trace(I^{-1})')
title('A-optimal objective function')
print(h, '-dpdf', 'fig2c.pdf')

% compute optimal constant flip angles 
[~, index_fD] = max(fD_vals);
[~, index_fE] = max(fE_vals);
[~, index_fA] = max(fA_vals);

fD_opt = 180/pi*theta_vals(index_fD)
fE_opt = 180/pi*theta_vals(index_fE)
fA_opt = 180/pi*theta_vals(index_fA)


%% Figure 3: Optimal flip angles for 2-site, 1-input model 

% specify optimization start point and options for MATLAB optimization toolbox 
initial_thetas_value = pi/2*ones(model.N, model.n + model.m);
options = optimset('MaxFunEvals', 20000, 'MaxIter', 500, 'Display', 'iter');

% perform optimization 
thetas_variable = optimal_flip_angle_design(model, design_criterion, ...
    initial_thetas_value, options); 

% plot optimal flip angles 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(thetas_variable.*180./pi, '.-', 'linewidth', 2, 'markersize', 20)  
title('Optimal flip angle scheme') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
legend('Pyr', 'Lac', 'AIF')
axis([1 model.N 0 100])
print(h, '-dpdf', 'fig3.pdf')


%% Figure 4: Comparison of flip angle schemes 

% compute D-optimal constant flip angle scheme 
options = optimset('MaxFunEvals', 20000, 'MaxIter', 500, 'Display', 'iter');
theta_constant = constant_optimal_flip_angle_design(model, ... 
    design_criterion, pi/2, options); 

% set flip angle values for the two constant flip angle schemes 
thetas_10 = 10/180*pi*ones(size(thetas_variable)); 
thetas_D = theta_constant*ones(size(thetas_variable)); 

% set number of observations of the random variable Y to sample 
num_runs = 20; 

% perform maximum likelihood estimation for each sample 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
hold on 
goodness_of_fit_criterion = 'maximum-likelihood'; 
j_max = ceil(length(model.parameters_of_interest)/2); 
for i = 1:num_runs
    i 
    
    % theta = 10 degrees
    y = generate_data(model, thetas_10); 
    [parameters_of_interest_est, nuisance_parameters_est] ...
       = parameter_estimation(y, model, goodness_of_fit_criterion, thetas_10); 
    save_10(i, :) = parameters_of_interest_est; 
    
    % constant D-optimal scheme
    y = generate_data(model, thetas_D); 
    [parameters_of_interest_est, nuisance_parameters_est] ...
       = parameter_estimation(y, model, goodness_of_fit_criterion, thetas_D); 
    save_D(i, :) = parameters_of_interest_est; 
    
    % variable D-optimal scheme 
    y = generate_data(model, thetas_variable); 
    [parameters_of_interest_est, nuisance_parameters_est] ...
       = parameter_estimation(y, model, goodness_of_fit_criterion, thetas_variable); 
    save_variable(i, :) = parameters_of_interest_est; 

end

% create scatterplot of the resulting parameter estimates 
h = figure;  
hold on 
plot(save_10(:, 1), save_10(:, 2), '.','markersize',40, 'Color', berkeley_colors(1,:))
plot(save_D(:, 1), save_D(:, 2), '.','markersize',40, 'Color', berkeley_colors(2,:))
plot(save_variable(:, 1), save_variable(:, 2), '.','markersize',40, 'Color', berkeley_colors(3,:))
plot(model.parameters_of_interest_nominal_values(1), model.parameters_of_interest_nominal_values(2), 'rx', 'markersize',20, 'linewidth',4)
legend('constant \theta = 10^\circ', ['D-optimal constant \theta = ', num2str(180/pi*theta_constant), '^\circ'], 'D-optimal variable', 'true model parameters' )
xlabel('k_{PL}')
ylabel('k_{TRANS}')
title('Spread of maximum-likelihood estimates')
hold off 
print(h, '-dpdf', 'fig4a.pdf')

% create bar graph of corresponding estimation errors 
h_fig = figure;  
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
h_bar_graph = bar( [mean(abs(save_10 - ones(num_runs, 1)*model.parameters_of_interest_nominal_values)); 
   mean(abs(save_D - ones(num_runs, 1)*model.parameters_of_interest_nominal_values)); 
   mean(abs(save_variable - ones(num_runs, 1)*model.parameters_of_interest_nominal_values))]'); 
for i = 1:length(h_bar_graph)
    set(h_bar_graph(i), 'FaceColor', berkeley_colors(i,:)) 
end
Labels = {'kPL', 'kTRANS'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels );
ylabel('mean error (1/s)') 
legend('constant \theta = 10^\circ', ['D-optimal constant \theta = ', num2str(180/pi*theta_constant), '^\circ'], 'D-optimal variable')
title('Mean error of maximum-likelihood estimates')
print(h_fig, '-dpdf', 'fig4b.pdf')


%% Figure 5: Initialization from different starting locations to demonstrate non-convexity 

for i=1:3
    i 
    % specify optimization start point and options for MATLAB optimization toolbox 
    if i == 1
        % for the first time, use 90 degrees as initialization 
        initial_thetas_value = pi/2*ones(model.N, model.n + model.m);
    else
        % after that, choose randomly 
        initial_thetas_value = pi/2*rand(model.N, model.n + model.m);
    end
    options = optimset('MaxFunEvals', 5000, 'MaxIter', 200, 'Display', 'iter');

    % perform optimization for 5000 steps 
    [thetas_variable, obj_val] = optimal_flip_angle_design(model, design_criterion, ...
        initial_thetas_value, options); 

    % plot optimal flip angles 
    h = figure; 
    set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
    plot(thetas_variable.*180./pi, '.-', 'linewidth', 2, 'markersize', 20)  
    title(['Optimal flip angle scheme -- objective value:' num2str(obj_val)]) 
    xlabel('acquisition number')
    ylabel('flip angle (degrees)')
    legend('Pyr', 'Lac', 'AIF')
    axis([1 model.N 0 100])
    if i == 1
        print(h, '-dpdf', 'fig5a.pdf')
    elseif i == 2
        print(h, '-dpdf', 'fig5c.pdf')
    elseif i == 3
        print(h, '-dpdf', 'fig5e.pdf')
    else
        error('should not be reachable')
    end
    
    
    options = optimset('MaxFunEvals', 20000, 'MaxIter', 500, 'Display', 'iter');

    % perform optimization for 20000 more steps 
    [thetas_variable, obj_val] = optimal_flip_angle_design(model, design_criterion, ...
        thetas_variable, options); 
    
    % plot optimal flip angles 
    h = figure;  
    set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
    plot(thetas_variable.*180./pi, '.-', 'linewidth', 2, 'markersize', 20)  
    title(['Optimal flip angle scheme -- objective value:' num2str(obj_val)]) 
    xlabel('acquisition number')
    ylabel('flip angle (degrees)')
    legend('Pyr', 'Lac', 'AIF')
    axis([1 model.N 0 100])
    if i == 1
        print(h, '-dpdf', 'fig5b.pdf')
    elseif i == 2
        print(h, '-dpdf', 'fig5d.pdf')
    elseif i == 3
        print(h, '-dpdf', 'fig5f.pdf')
    else
        error('should not be reachable')
    end
    
end 

%% Figure 7: Optimal flip angles for 2-site model with constraint that pyruvate flip angles applied to input and first compartment are equal 

% set constraints on input 
model.flip_angle_input_matrix = [1 0; 
                                 0 1; 
                                 1 0]; 

% specify optimization start point and options for MATLAB optimization toolbox 
initial_thetas_value = pi/2*ones(model.N, model.n);
options = optimset('MaxFunEvals', 20000, 'MaxIter', 500, 'Display', 'iter');

% perform optimization 
thetas_variable = optimal_flip_angle_design(model, design_criterion, ...
    initial_thetas_value, options); 

% plot optimal flip angles 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(thetas_variable.*180./pi, '.-', 'linewidth', 2, 'markersize', 20)  
title('Optimal flip angle scheme') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
legend('Pyr', 'Lac', 'AIF')
axis([1 model.N 0 100])
print(h, '-dpdf', 'fig7.pdf')


%% Figure 6: Optimal flip angles for 3-site exchange model 

clear model 

% initialize model object 
model = linear_exchange_model; 

% define model parameters
syms R1P R1L R1A kPL kPA kTRANS 
% define input parameters 
syms t0 alpha_1 beta_1 A0 
% define noise parameters 
syms sigma_1 sigma_2 sigma_3
% define initial state parameters
syms P0 L0 

% parameters of interest 
% (those for which we wish to compute an estimate with minimal variance) 
model.parameters_of_interest = [kPL kPA kTRANS]; 
model.parameters_of_interest_nominal_values = [0.02 0.03 0.04]; 

% nuisance parameters
% (those parameters that are unknown but whose estimates we only care about
% insofar as they allow us to estamate the parameters of interest) 
model.nuisance_parameters = [alpha_1 beta_1 A0];
model.nuisance_parameters_nominal_values = [ 2  5  1]; 

% known parameters
% (those whose values are assumed to be known constants) 
model.known_parameters = [R1P R1L R1A t0 P0 L0]; 
model.known_parameter_values = [1/35 1/30 1/40 0 0 0];

% three-site exhange model 
model.A = [ -kPL-kPA-R1P  0    0;
            kPL         -R1L  0;
            kPA          0   -R1A];   
model.B = [kTRANS; 0; 0]; 

% define input function shape  
model.u = @(t) A0 * (t - t0)^alpha_1 *exp(-(t - t0)/beta_1); % gamma-variate input  

% define initial condition 
model.x0 = [P0; L0; 0]; 

% define repetition time
model.TR = 2; 

% define number of acquisitions 
model.N = 25; 

% choose noise type
model.noise_type = 'Rician';

% choose noise magnitude  
model.noise_parameters = [0.01 0.01 0.01 0.1]; 

% choose flip angle input matrix 
model.flip_angle_input_matrix = eye(model.m + model.n);                               

% choose design criterion 
design_criterion = 'D-optimal'; 

% discretize model (doing this in advance makes things run faster) 
model = discretize(model);  

% compute sensitivities (doing this in advance makes things run faster)
if ~model.sensitivities_computed ...
        && (strcmp(design_criterion, 'D-optimal') ...
            || strcmp(design_criterion, 'E-optimal') ...
            || strcmp(design_criterion, 'A-optimal') ...
            || strcmp(design_criterion, 'T-optimal') ...
           )
    model = sensitivities(model);  
end


% specify optimization start point and options for MATLAB optimization toolbox 
initial_thetas_value = pi/2*ones(model.N, model.n + model.m);
options = optimset('MaxFunEvals', 20000, 'MaxIter', 500, 'Display', 'iter');

% perform optimization 
thetas_variable = optimal_flip_angle_design(model, design_criterion, ...
    initial_thetas_value, options); 

% plot optimal flip angles 
h = figure; 
set(gca,'ColorOrder', berkeley_colors, 'NextPlot', 'replacechildren')
plot(thetas_variable.*180./pi, '.-', 'linewidth', 2, 'markersize', 20)  
title('Optimal flip angle scheme') 
xlabel('acquisition number')
ylabel('flip angle (degrees)')
legend('Pyr', 'Lac', 'Ala', 'AIF')
axis([1 model.N 0 100])
print(h, '-dpdf', 'fig6.pdf')
