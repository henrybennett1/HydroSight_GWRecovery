classdef derivedweighting_PearsonsPositiveRescaled < derivedResponseFunction_abstract
% Pearson's type III impulse response transfer function class. 

    properties(GetAccess=public, SetAccess=protected)
        A
        settings 
    end
%%  STATIC METHODS        
% Static methods used by the Graphical User Interface to inform the
% user of the available model options and their input format.
    methods(Static)        
        function [modelSettings, colNames, colFormats, colEdits,tooltipString] = modelOptions(bore_ID, forcingDataSiteID, siteCoordinates)
           modelSettings = {};
           colNames = {};
           colFormats = {};
           colEdits = [];       
           tooltipString='';
        end
        function modelDescription = modelDescription()
           modelDescription = {'Name: derivedweighting_PearsonsPositiveRescaled', ...
                               '', ...
                               'Purpose: simulation of recharge-like climate forcing (ie inputs of rainfall, free drainage etc) using previously created', ...
                               'responseFunction_Pearsons weighting function. ', ...
                               '', ...                               
                               'Number of parameters: 1', ...
                               '', ...                               
                               'Options: none', ...
                               '', ...                               
                               'Comments: In combination with the derived forcing function derivedForcing_linearUnconstrainedScaling, this function can be used ', ...
                               'to estimate the impact of revegetation on recharge by weighting a recharge output from the aformentioned function (which needs an', ...
                               'input for the landuse change).', ...
                               '', ...                               
                               'References: (none)'};
        end        
    end
%%    
    methods
        % Constructor
        function obj = derivedweighting_PearsonsPositiveRescaled(bore_ID, forcingDataSiteID, siteCoordinates, sourceResponseFunctionObject, options, params)
                        
            % Define default parameters 
            if nargin==5
                params=0.01;
            end
                
            % Set parameters for transfer function.
            setParameters(obj, params)     
            
            % Assign the source object to settings.
            obj.settings.sourceObject = sourceResponseFunctionObject;
            
            % Check the source model can return a normalised theta value. 
            % The normalised theta results is rescaled within this function.
            % This is undertaken to reduce parameter covariance.            
            if ~any(strcmp(methods(obj.settings.sourceObject),'theta_normalised'))
                error('This weighting function normalises and rescales the Pearsons weigthing function. However, the following expected method within the source function could not be found: "theta_normalised"');
            end

        end
       
        % Set parameters
        function setParameters(obj, params)
            obj.A = params(1,:);
        end
        
        % Get model parameters
        function [params, param_names] = getParameters(obj)
            params(1,:) = obj.A;
            param_names = {'A'};        
        end        
        
        function isValidParameter = getParameterValidity(obj, params, param_names)                                    
            % Get physical bounds.
            [params_upperLimit, params_lowerLimit] = getParameters_physicalLimit(obj);

    	    % Check parameters are within bounds.
            isValidParameter = params >= params_lowerLimit(:,ones(1,size(params,2))) & ...
                    params <= params_upperLimit(:,ones(1,size(params,2)));
        end
        
        % Return fixed upper and lower bounds to the parameters.
        function [params_upperLimit, params_lowerLimit] = getParameters_physicalLimit(obj)
            params_upperLimit = log10(1.0);
            params_lowerLimit = log10(eps());

        end        
        
        % Return fixed upper and lower plausible parameter ranges. 
        % This is used to define reasonable range for the initial parameter sets
        % for the calibration. These parameter ranges are only used in the 
        % calibration if the user does not input parameter ranges.
        function [params_upperLimit, params_lowerLimit] = getParameters_plausibleLimit(obj)
            params_upperLimit = log10(0.1);
            params_lowerLimit = log10(sqrt(eps()));
        end
        
        % Calculate impulse-response function.
        function result = theta(obj, t)           
            % Call the source model theta function
            result = theta_normalised(obj.settings.sourceObject, t);
            
            % Scale the result by back transformed 'A'.
            result = 10.^obj.A .* result;            
        end   
        
        % Calculate integral of impulse-response function from t to inf.
        % This is used to minimise the impact from a finit forcign data
        % set.
        function result = intTheta_upperTail2Inf(obj, t)                       
            % Get the original A parameter scaling value.
            [result, A_orig] = theta_normalised(obj.settings.sourceObject, 1);
            
            % Call the source model theta function
            result = intTheta_upperTail2Inf(obj.settings.sourceObject, t);            
            
            % Rescale the result.
            result = 10.^obj.A .* (result./A_orig);
        end   

        % Calculate integral of impulse-response function from 0 to 1.
        % This is used handle rapidly chnageing fucntion in the range from 0 to 1.        
        function result = intTheta_lowerTail(obj, t)  
            % Get the original A parameter scaling value.
            [result, A_orig] = theta_normalised(obj.settings.sourceObject, 1);            
            
            % Call the source model theta function
            result = intTheta_lowerTail(obj.settings.sourceObject, t);
            
            % Rescale the result.
            result = 10.^obj.A .* (result./(10.^obj.A));
        end
                
        % Transform the estimate of the response function * the forcing.
        function result = transform_h_star(obj, h_star_est)
           result = h_star_est(:,end);
        end   
        
        % Return the derived variables.
        function [params, param_names] = getDerivedParameters(obj)
            params(1,:) = 10.^obj.A;
            param_names = {'A : Pearsons weighting function rescaling'};                
        end

        function delete(obj)
% delete class destructor
%
% Syntax:
%   delete(obj)
%
% Description:
%   Loops through parameters and, if not an object, empties them. Else, calls
%   the sub-object's destructor.
%
% Input:
%   obj -  model object
%
% Output:  
%   (none)
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure Engineering, 
%   The University of Melbourne.
%
% Date:
%   24 Aug 2016
%%            
            propNames = properties(obj);
            for i=1:length(propNames)
               if isempty(obj.(propNames{i}))
                   continue;
               end                
               if isobject(obj.(propNames{i}))
                delete(obj.(propNames{i}));
               else               
                obj.(propNames{i}) = []; 
               end
            end
        end            
        
    end

end

