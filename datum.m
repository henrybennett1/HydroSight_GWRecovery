classdef datum
    %DATUM Summary of this class goes here
    %   Detailed explanation goes here
    %   datum_shift is the distance from the minimum head value that we
    %   want to define the datum to be, this can be defined by a range of
    %   methods that are TBD.

    %   model_shift is the difference between the default model and the
    %   secondary model that is being made to compare the two states
    
    properties
        d
        params_upperLimit
        params_lowerLimit
    end
    
    methods
        function obj = datum(obs_head)
            obj.d = min(obs_head(:,end));
            obj.params_upperLimit = obj.d + 50;
            obj.params_lowerLimit = obj.d - 50;
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
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
        end
        
        function [params_upperLimit, params_lowerLimit] =  getParameters_physicalLimit(obj)      
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
        end    
    end
end

