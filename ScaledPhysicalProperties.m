classdef ScaledPhysicalProperties
    properties
        nu
        kappa
        rho
        cp
        k
        alpha_g
        H_halflife
        H_0

        unscaled_properties
        scaling
        T_offset
    end
    methods
        function obj=ScaledPhysicalProperties(unscaled_properties, scaling)
            obj.nu = unscaled_properties.nu / scaling.nu;

            obj.kappa = unscaled_properties.kappa / scaling.kappa;
            obj.rho = unscaled_properties.rho / scaling.rho;
            obj.cp = unscaled_properties.cp / scaling.cp;
            obj.k = obj.kappa * obj.rho * obj.cp;

            obj.alpha_g = unscaled_properties.alpha_g / scaling.alpha_g;
            obj.H_halflife = unscaled_properties.H_halflife / scaling.t;
            obj.H_0 = unscaled_properties.H_0 / scaling.H;
            obj.unscaled_properties = unscaled_properties;
            obj.scaling = scaling;
        end
    end
end
