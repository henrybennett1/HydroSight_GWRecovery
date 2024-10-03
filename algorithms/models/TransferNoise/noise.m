classdef noise < handle 
    %NOISE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        alpha
        sigma_n
        params_upperLimit
        params_lowerLimit   %Do we need these values here?
        variables
        
    end
    
    methods
        function obj = noise(timesteps)
            obj.alpha = log10(0.1);
            obj.sigma_n = 0.1;

            delta_time = diff(timesteps);

            % This parameter is assumed to be the noise parameter 'alpha'.  
            alpha_upperLimit = 5; 
            % while abs(sum( exp( -2.*alpha_upperLimit .* delta_time ) )) < eps() ...
            % || exp(mean(log( 1- exp( -2.*alpha_upperLimit .* delta_time) ))) < eps()
            %     alpha_upperLimit = alpha_upperLimit - 0.01;
            %     if alpha_upperLimit <= eps()                                   
            %         break;
            %     end
            % end
            % if alpha_upperLimit <= eps()
            %     alpha_upperLimit = inf;
            % else
            %     % Transform alpha log10 space.
            %     alpha_upperLimit = log10(alpha_upperLimit);
            % end                           
                        

            obj.params_upperLimit =  alpha_upperLimit;
            % obj.params_lowerLimit = log10(sqrt(eps()))+4;
            obj.params_lowerLimit = 0;
            obj.variables.delta_time = delta_time;
            
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
        function noise = getNoise(obj, time_points, noisePercnile)
            
            % Set percentile for noise
            if nargin==2
                noisePercnile = 0.95;
            else
                noisePercnile = noisePercnile(1);
            end
            
            noise = obj.sigma_n; %repath for each sigma values

            nparamsets = length(noise);
            if nparamsets>1
                noise = reshape(noise, 1, 1, nparamsets);
            end
            
            noise = [repmat(time_points,1,1,nparamsets), ones(size(time_points,1),2) .* norminv(noisePercnile,0,1).* noise]; 
            %Why column size of 2?
        end
    end
end

