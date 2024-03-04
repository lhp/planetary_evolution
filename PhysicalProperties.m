classdef PhysicalProperties
    properties
        nu
        kappa
        rho
        cp
        k
        alpha_g
        H_halflife
        H_0  % Heating in W/kg

    end
    methods
        function obj=PhysicalProperties(nu, rho, cp, k, alpha_g, H_halflife, H_0)
            obj.nu = nu;
            
            obj.rho = rho;
            obj.cp = cp;
            obj.k = k;
            obj.kappa = k / (rho*cp);

            obj.alpha_g = alpha_g;
            obj.H_halflife = H_halflife;
            obj.H_0 = H_0;
        end
    end
end
