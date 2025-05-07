% E.J. VanHenten andG. VanStraten, Sensitivity analysis of a dynamic growth model of lettuce, Journal of Agricultural Engineering Research, vol. 59(1), pp. 19–31, 1994.
function params = load_parameters()
    params.c_alpha      = 0.68;
    params.c_beta       = 0.8;
    params.c_gr_max     = 5.0e-6;
    params.c_gamma      = 1.0;
    params.c_Q10_gr     = 1.6;
    params.c_resp_sht   = 3.47e-7;
    params.c_resp_rt    = 1.16e-7;
    params.c_tau        = 0.07;
    params.c_Q10_resp   = 2.0;
    params.c_K          = 0.9;
    params.c_lar        = 75.0e-3;
    params.c_epsilon    = 17.0e-6;
    params.c_Gamma      = 40;
    params.c_10_Gamma   = 2.0;
    params.c_bnd        = 0.004;
    params.c_stm        = 0.02;
    params.c_car1       = -1.32e-5;
    params.c_car2       = 5.94e-4;
    params.c_car3       = -2.64e-3;
    params.c_w          = 7.32e-3;
end