classdef datum < handle
    %DATUM Summary of this class goes here
    %   Detailed explanation goes here
    %   datum_shift is the distance from the minimum head value that we
    %   want to define the datum to be, this can be defined by a range of
    %   methods that are TBD.

    %   model_shift is the difference between the default model and the
    %   secondary model that is being made to compare the two states
    
    properties
        d
        variables
    end
    
    methods
        function obj = datum(obs_head)
            %obj.d = min(obs_head(:,end)) - 0.1 * abs(min(obs_head(:,end)));
            obj.d = min(obs_head(:,end));
            obj.variables.upper_phys_d = min(obs_head(:,end))*1.2;
            obj.variables.lower_phys_d = min(obs_head(:,end))*0.8;
            obj.variables.upper_plaus_d = min(obs_head(:,end))*1.1;
            obj.variables.lower_plaus_d = min(obs_head(:,end))*0.9;
            % Potentially include a modifier term that changes the datum so
            % they're different from each other
            %Datum shift should be defined based on the drainage elevation
            %value that is found in model_TFN line 1868
        end
        
        function d = getdatum(obj)
            d = obj.d;
        end

        function [params, param_names] = getParameters(obj)
            params = getdatum(obj);
            param_names = {'d'};
        end
        
        function [params_upperLimit, params_lowerLimit] =  getParameters_plausibleLimit(obj)      
            params_upperLimit = obj.variables.upper_plaus_d;
            params_lowerLimit = obj.variables.lower_plaus_d;
        end
        
        function [params_upperLimit, params_lowerLimit] =  getParameters_physicalLimit(obj)      
            params_upperLimit = obj.variables.upper_phys_d;
            params_lowerLimit = obj.variables.lower_phys_d;
        end 

        function setParameters(obj, params)
            obj.d = params(1,:);   
        end

        function isValidParameter = getParameterValidity(obj, params, param_names)
            [params_upperLimit, params_lowerLimit] = getParameters_physicalLimit(obj);

            % Check parameters are within bounds.
            isValidParameter = params >= params_lowerLimit(:,ones(1,size(params,2))) & ...
    		params <= params_upperLimit(:,ones(1,size(params,2)));
        end


    end
end

