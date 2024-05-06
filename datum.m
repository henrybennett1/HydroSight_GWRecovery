classdef datum
    %DATUM Summary of this class goes here
    %   Detailed explanation goes here
    %   datum_shift is the distance from the minimum head value that we
    %   want to define the datum to be, this can be defined by a range of
    %   methods that are TBD.

    %   model_shift is the difference between the default model and the
    %   secondary model that is being made to compare the two states
    
    properties
        d1
        d2
        params_upperLimit
        params_lowerLimit
    end
    
    methods
        function obj = datum(obs_head,model_shift,datum_shift)
            obj.d1 = min(obs_head) - datum_shift;                               %DEFAULT MODEL
            obj.d2 = min(obs_head) - datum_shift - model_shift;                 %SHIFTED MODEL
            obj.params_upperLimit = 1000;
            obj.params_lowerLimit = -1000;
        end
        
        function [d1,d2] = getdatum(obj)
            d1 = obj.d1;
            d2 = obj.d2;
        end
        
        function [params_upperLimit, params_lowerLimit] =  getParameters_plausibleLimit(obj)      
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
        end
    end
end

