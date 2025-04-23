function dXdt = lettuceODE(t, X, params, inputs)
% Lettuce growth model.
%
% Args:
%   t: Current time (not explicitly used here as inputs are constant, but required by ODE solvers)
%   X: State vector [X_nsdw; X_sdw] (g m^-2)
%   params: Structure containing model parameters
%   inputs: Structure containing environmental inputs
%
% Returns:
%   dXdt: Column vector of derivatives [dX_nsdw/dt; dX_sdw/dt]

% Unpack states
X_nsdw = X(1); % Non-structural dry weight (g m^-2)
X_sdw = X(2);  % Structural dry weight (g m^-2)

% Ensure non-negative weight
X_nsdw = max(0, X_nsdw);
X_sdw = max(eps, X_sdw); % Avoid division by zero

% Parameters
c_alpha      = params.c_alpha;     % CO2 to CH2O conversion factor
c_beta       = params.c_beta;      % Growth efficiency factor
c_gr_max     = params.c_gr_max;    % Max specific growth rate (s^-1)
c_gamma      = params.c_gamma;     % Growth rate balance factor
c_Q10_gr     = params.c_Q10_gr;    % Q10 factor for growth
c_resp_sht   = params.c_resp_sht;  % Shoot maint. respiration coeff (s^-1) at 25C
c_resp_rt    = params.c_resp_rt;   % Root maint. respiration coeff (s^-1) at 25C
c_tau        = params.c_tau;       % Root dry weight fraction
c_Q10_resp   = params.c_Q10_resp;  % Q10 factor for respiration
c_K          = params.c_K;         % Canopy light extinction coeff
c_lar        = params.c_lar;       % Structural leaf area ratio (m^2 g^-1)
c_epsilon   = params.c_epsilon;  % Light use efficiency coeff (g J^-1) at high CO2
c_Gamma     = params.c_Gamma;    % CO2 compensation point (ppm) at 20C
c_10_Gamma   = params.c_10_Gamma;  % Q10 factor for Gamma
c_bnd        = params.c_bnd;       % Boundary layer conductance (m s^-1)
c_stm        = params.c_stm;       % Stomatal conductance (m s^-1)
c_car1       = params.c_car1;      % Carboxylation conductance param 1
c_car2       = params.c_car2;      % Carboxylation conductance param 2
c_car3       = params.c_car3;      % Carboxylation conductance param 3
c_w          = params.c_w;         % Density of CO2 (g m^-3)

U_PAR = inputs.U_PAR; % Incident PAR (W m^-2)
U_CO2 = inputs.U_CO2; % CO2 concentration (ppm)
U_T   = inputs.U_T;   % Canopy temperature (degrees C)

% Intermediate calculations


% Light use efficiency (epsilon)
Gamma = c_Gamma * c_10_Gamma^((U_T - 20) / 10);
if U_CO2 + 2*Gamma <= 0 % prevents division by zero
    epsilon = 0;
else
    epsilon = c_epsilon * (U_CO2 - Gamma) / (U_CO2 + 2*Gamma); 
    epsilon = max(0, epsilon); % cannot be negative
end

% Carboxylation conductance (1/g_CO2)
c_car = c_car1 * U_T^2 + c_car2 * U_T + c_car3; 
c_car = max(eps, c_car); 

g_CO2_inv = (1 / max(eps, c_bnd)) + (1 / max(eps, c_stm)) + (1 / max(eps, c_car));
g_CO2 = 1 / g_CO2_inv; 

% Max Canopy Photosynthesis (f_phot_max)
phot_denom = epsilon * U_PAR + g_CO2 * c_w * (U_CO2 - Gamma);
if phot_denom <= 0
    f_phot_max = 0;
else
    f_phot_max = (epsilon * U_PAR * g_CO2 * c_w * (U_CO2 - Gamma)) / phot_denom;
    f_phot_max = max(0, f_phot_max); % cannot be negative
end

% Gross Canopy Photosynthesis Rate (f_phot)
f_phot = (1 - exp(-c_K * c_lar * (1 - c_tau) * X_sdw)) * f_phot_max;

% Growth Rate (r_gr)
growth_denom = c_gamma * X_sdw + X_nsdw;
if growth_denom <= 0 || X_nsdw < 0
    r_gr = 0;
else
    r_gr = c_gr_max * (X_nsdw / growth_denom) * c_Q10_gr^((U_T - 20) / 10); 
    r_gr = max(0, r_gr); % cannot be negative
end

% Maintenance Respiration (f_resp)
f_resp = (c_resp_sht * (1 - c_tau) * X_sdw + c_resp_rt * c_tau * X_sdw) * c_Q10_resp^((U_T - 25) / 10); 
f_resp = max(0, f_resp); % cannot be negative

% dX_nsdw/dt
dX_nsdw_dt = c_alpha * f_phot - r_gr * X_sdw - f_resp - ((1 - c_beta) / c_beta) * r_gr * X_sdw;

% dX_sdw/dt
dX_sdw_dt = r_gr * X_sdw;

% Return derivatives as a col vector
dXdt = [dX_nsdw_dt; dX_sdw_dt];

end