function [model, scaling, total_ticks]=setup_model(...
    d, H, normalized_core_C, Ra, L_l, kappa_l, solid, total_time, ...
    velocity_skip, settling_tolerance)

    gy_in_s = 1e9*365*24*60*60;
    L_r = sum(d);
    T_surface = 400;
    dT_r = 1500;
    alphag_r = 2.5e-5*3.7;
    rho_r = 3200;
    cp_r = 1200;
    k_r = 4;
    kappa_r = k_r / (rho_r * cp_r);
    nu_r = alphag_r * dT_r * L_r^3 / (kappa_r * Ra);
    H_halflife_r = 1.75 * gy_in_s;
    pp_r = PhysicalProperties(nu_r, rho_r, cp_r, k_r, alphag_r, H_halflife_r, H);
    
    dT_l = 1;
    nu_l = 1;
    rho_l = 1;
    cp_l = 1;
    scaling = Scaling(L_r/L_l, kappa_r/kappa_l, nu_r/nu_l, rho_r/rho_l, cp_r/cp_l, dT_r/dT_l, T_surface);
    pp_l = ScaledPhysicalProperties(pp_r, scaling);
    
    ny = L_l;
    ny_convective = L_l * max(d)/sum(d);
    nx = 4 * round(sqrt(2) * ny_convective);
    
    
    dTdt_r_0 = H / cp_r;  % W/kg=J/kg/s / J/kg/K = K/s
    [solid_grid, dTdt_l_0_grid] = fn_get_grid_properties_from_layers(...
        nx, ny, d/scaling.x, solid, pp_l.H_0 / pp_l.cp);
    
    T_0_l = fn_create_init_T(nx, ny, dT_l, 0);
    core_T_0_l = dT_l;  % No overheating.
    
    total_ticks = round(total_time / scaling.t);
    
    model = PlanetModel(T_0_l, core_T_0_l, normalized_core_C, ...
                    pp_l, dTdt_l_0_grid, solid_grid, total_ticks, velocity_skip, settling_tolerance);
end


%% functions

function [solid, dTdt]=fn_get_grid_properties_from_layers(...
    nx, ny, layer_d, layer_solid, layer_dTdt)

    solid = zeros(nx, ny);
    dTdt = zeros(nx, ny);

    layer_boundaries = round([0 cumsum(layer_d)]);
    % Check they are all integers.
    for i=1:length(layer_boundaries)
        if layer_boundaries(i) ~= int64(layer_boundaries(i))
            fprintf("Error: layer %g does not fall on a lattice point.\n", i);
        end
    end
    
    for i=1:length(layer_boundaries)-1
        solid(:, layer_boundaries(i)+1:layer_boundaries(i+1)) = layer_solid(i);
        dTdt(:, layer_boundaries(i)+1:layer_boundaries(i+1)) = layer_dTdt(i);
    end

    % Do this last to ensure that we have solid top no matter what.
    solid(:, end, 1) = 1;
end

function T=fn_create_init_T(nx, ny, T_bot, T_top)
    T = nan(nx, ny, 1);
    T_background = nan(size(T));
    lambda = round(nx/(sqrt(2)*ny));
    perturbation = (T_bot-T_top)/10000 * sin(2*pi*lambda*(1:nx)/nx);
    for j=1:ny
        T_background(:, j, 1) = T_bot;  % Initial condition: all hot.
        T(:, j, 1) = T_background(:, j, 1) + perturbation';
    end
end
