function model = lettuce_growth_model()
% Creates the FULL symbolic representation of the lettuce growth model.

    % ------- Define symbolic variables -------

    % State Variables
    syms X_nsdw X_sdw

    % Parameters
    syms c_alpha c_beta c_gr_max c_gamma c_Q10_gr ...        % Growth/Efficiency 
         c_resp_sht c_resp_rt c_tau c_Q10_resp ...           % Respiration 
         c_K c_lar ...                                       % Light interception 
         c_epsilon_max c_Gamma_base c_10_Gamma ...           % Photosynthesis efficiency
         c_bnd c_stm c_car1 c_car2 c_car3 ...                % CO2 conductance components 
         c_w                                                 % CO2 density 

    % Inputs
    syms U_PAR U_CO2 U_T

    % Initial Conditions
    syms X_nsdw0 X_sdw0

    % ------- Intermediate calculations (using full symbolic parameters) -------

    % Gamma
    Gamma = c_Gamma_base * c_10_Gamma^((U_T - 20)/10);
    CO2_driving_force = U_CO2 - Gamma;

    % epsilon
    epsilon_den = U_CO2 + 2*Gamma;
    epsilon = c_epsilon_max * CO2_driving_force / epsilon_den; 

    c_car = c_car1 * U_T^2 + c_car2 * U_T + c_car3;

    g_CO2_inv = (1 / c_bnd) + (1 / c_stm) + (1 / c_car);
    g_CO2 = 1 / g_CO2_inv;

    % 5f_phot_max
    phot_term_a = epsilon * U_PAR;
    phot_term_b = g_CO2 * c_w * CO2_driving_force;
    phot_max_den = phot_term_a + phot_term_b;
    f_phot_max = (phot_term_a * phot_term_b) / phot_max_den; 

    % f_phot
    light_interception = 1 - exp(-c_K * c_lar * (1 - c_tau) * X_sdw);
    f_phot = light_interception * f_phot_max;

    % r_gr
    growth_rate_ratio_den = c_gamma * X_sdw + X_nsdw;
   growth_rate_ratio = X_nsdw / growth_rate_ratio_den; 
    r_gr = c_gr_max * growth_rate_ratio * c_Q10_gr^((U_T - 20)/10);

    % f_resp
    f_resp = (c_resp_sht*(1 - c_tau)*X_sdw + c_resp_rt*c_tau*X_sdw) * c_Q10_resp^((U_T - 25)/10); % Uses symbolic base rates

    % ------- Model Structure -------

    % States
    model.sym.x = [X_nsdw; X_sdw];

    % Parameters
    model.sym.p = [c_alpha; c_beta; c_gr_max; c_gamma; c_Q10_gr; ...
                   c_resp_sht; c_resp_rt; c_tau; c_Q10_resp; ...
                   c_K; c_lar; ...
                   c_epsilon_max; c_Gamma_base; c_10_Gamma; ...
                   c_bnd; c_stm; c_car1; c_car2; c_car3; ...
                   c_w];

    % Inputs
    model.sym.u = [U_PAR; U_CO2; U_T];

    % Initial conditions
    model.sym.x0 = [X_nsdw0; X_sdw0];

    % Outputs
    model.sym.y = [X_nsdw; X_sdw];
    model.sym.g = [0; 0];

    % Right-hand side function 
    growth_term = r_gr * X_sdw;
    growth_loss_term = ((1 - c_beta)/c_beta) * growth_term; % Assumes 0 < c_beta < 1
    dX_nsdw_dt = c_alpha * f_phot - growth_term - f_resp - growth_loss_term;
    dX_sdw_dt  = growth_term;

    model.sym.xdot = [dX_nsdw_dt; dX_sdw_dt];

    model.param_names = {'c_alpha', 'c_beta', 'c_gr_max', 'c_gamma', 'c_Q10_gr', ...
                         'c_resp_sht', 'c_resp_rt', 'c_tau', 'c_Q10_resp', ...
                         'c_K', 'c_lar', ...
                         'c_epsilon_max', 'c_Gamma_base', 'c_10_Gamma', ...
                         'c_bnd', 'c_stm', 'c_car1', 'c_car2', 'c_car3', ...
                         'c_w'};
    model.state_names = {'X_nsdw', 'X_sdw'};
    model.input_names = {'U_PAR', 'U_CO2', 'U_T'};
    model.output_names = {'y_X_nsdw', 'y_X_sdw'};

end