classdef probabilities
    %PROBABILITIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        initial
        trans_state1
        trans_state2
        params_upperLimit
        params_lowerLimit
    end
    
    methods
        function obj = probabilities(initial,trans1,trans2)
            obj.initial = initial;
            obj.trans_state1 = trans1;
            obj.trans_state2 = trans2;
            obj.params_upperLimit = 1;
            obj.params_lowerLimit = 0;
        end
        
        function [initial1, initial2] = getinitial_probabilities(obj)
            initial1 = obj.initial;
            initial2 = 1 - obj.initial;
        end

        function [trans, maintain] = gettrans_probabilities(obj,state)
            if state == 1
                trans = obj.trans_state1;
                maintain = 1 - obj.trans_state1;
            elseif state == 2
                trans = obj.trans_state2;
                maintain = 1 - obj.trans_state2;
            else
                error("State value not properly defined within bounds");
            end
        end
        function [params_upperLimit, params_lowerLimit] =  getParameters(obj)      
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
            %FIX
        end
        function [params_upperLimit, params_lowerLimit] =  getParameters_plausibleLimit(obj)
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
        end
        function [params_upperLimit, params_lowerLimit] =  getParameters_physicalLimit(obj)
            params_upperLimit = obj.params_upperLimit;
            params_lowerLimit = obj.params_lowerLimit;
            %FIX
        end
    end
end

