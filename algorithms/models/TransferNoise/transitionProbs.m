classdef transitionProbs < handle
    %PROBABILITIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        initial
        trans_state1
        trans_state2
    end
    
    methods
        function obj = transitionProbs(initial,trans1,trans2)
            obj.initial = initial;
            obj.trans_state1 = trans1;
            obj.trans_state2 = trans2;
        end
        function [params, param_names] = getParameters(obj)
            params = [obj.initial;obj.trans_state1;obj.trans_state2];
            param_names = {'Inital Prob';'Transition State 1'; 'Transition State 2'};

        end
        function setParameters(obj, params)
            obj.initial = params(1,:);
            obj.trans_state1 = params(2,:);
            obj.trans_state2 = params(3,:);
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
        function [params_upperLimit, params_lowerLimit] =  getParameters_plausibleLimit(obj)
            params_upperLimit = [1;1;1];
            params_lowerLimit = [0;0;0];
       
        end

        function [params_upperLimit, params_lowerLimit] =  getParameters_physicalLimit(obj)
            params_upperLimit = [1;1;1];
            params_lowerLimit = [0;0;0];
            
        end 
        function isValidParameter = getParameterValidity(obj, params, param_names)
            [params_upperLimit, params_lowerLimit] = getParameters_physicalLimit(obj);
            isValidParameter = params >= params_lowerLimit(:,ones(1,size(params,2))) & ...
                    params <= params_upperLimit(:,ones(1,size(params,2)));
        end
    end
end

