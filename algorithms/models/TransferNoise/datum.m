classdef datum < handle
    %DATUM Summary of this class goes here
    %   Detailed explanation goes here
    %   datum_shift is the distance from the minimum head value that we
    %   want to define the datum to be, this can be defined by a range of
    %   methods that are TBD.

    %   model_shift is the difference between the default model and the
    %   secondary model that is being made to compare the two states
    
    properties
        headDatum
        params_upperLimit
        params_lowerLimit
    end
    
    methods
        function obj = datum(obs_head)
            obj.headDatum = min(obs_head(:,end));                               %DEFAULT MODEL
            obj.params_lowerLimit = min(obs_head(:,end)) *0.5;
            obj.params_upperLimit = max(obs_head(:,end)) *1.5;

        end
        function [params, param_names] = getParameters(obj)
            params = obj.headDatum;
            param_names = {'Head Datum'};

        end

        function [d1] = getdatum(obj)
            d1 = obj.headDatum;
            
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

        %function [params_upperLimit, params_lowerLimit] =  getParameters_plausibleLimit(obj)      
            %DO WE NEED TO MAKE ONE FOR THESE?
        %end
    end
end

