classdef noise < handle 
    %NOISE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma_n
        params_upperLimit
        params_lowerLimit
    end
    
    methods
        function obj = noise(timesteps)
            %NOISE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sigma_n = 0.1;

            delta_time = diff(timesteps);

            % This parameter is assumed to be the noise parameter 'alpha'.  
            alpha_upperLimit = 100; 
            while abs(sum( exp( -2.*alpha_upperLimit .* delta_time ) )) < eps() ...
            || exp(mean(log( 1- exp( -2.*alpha_upperLimit .* delta_time) ))) < eps()
                alpha_upperLimit = alpha_upperLimit - 0.01;
                if alpha_upperLimit <= eps()                                   
                    break;
                end
            end
            if alpha_upperLimit <= eps()
                alpha_upperLimit = inf;
            else
                % Transform alpha log10 space.
                alpha_upperLimit = log10(alpha_upperLimit);
            end                           
                        

            obj.params_upperLimit =  alpha_upperLimit;
            obj.params_lowerLimit = log10(sqrt(eps()))+4;
            
        end
        
        function [params, param_names] = getParameters(obj)
            params = obj.sigma_n;
            param_names = {'sigma_n'};

        end
        function setParameters(obj, params)
            obj.sigma_n = params(1,:);   

        end
        function [params_upperLimit, params_lowerLimit] =  getParameters_plausibleLimit(obj)
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
        end

        function [params_upperLimit, params_lowerLimit] =  getParameters_physicalLimit(obj)
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
            %FIX and add to all others
        end  

        function isValidParameter = getParameterValidity(obj, params, param_names)
            [params_upperLimit, params_lowerLimit] = getParameters_physicalLimit(obj);
            isValidParameter = params >= params_lowerLimit(:,ones(1,size(params,2))) & ...
                    params <= params_upperLimit(:,ones(1,size(params,2)));
        end

    end
end

