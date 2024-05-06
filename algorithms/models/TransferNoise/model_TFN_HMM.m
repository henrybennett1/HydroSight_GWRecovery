classdef model_TFN_HMM < model_TFN
% Class definition for Transfer Function Noise (TFN) model for use with HydroSight
%
% Description: 
%   The class defines a transfer function Noise (TFN) model for
%   simulating time series of groundwater head. The model should be
%   defined, calibrated and solved using HydroSight() or 
%   the graphical user interface.
%
%   This class uses an object-oriented structure to provide a highly flexible 
%   model structure and the efficient inclusion of new structural
%   componants. Currently, linear and nonlinear TFN model can be built.
%   The linear models are based upon von Asmuth et al. 2002 and von Asmuth et al. 
%   2008 and the nonlinear models use a nonlinear transform of the input
%   climate data to account for runoff and nonlinear free drainage (see
%   Peterson & Western, 2014) 
%
%   Additional key features of the model include:
%
%       - long term historic daily climate forcing is used to estimate the
%       groundwater head at each time point via use of continuous transfer
%       function. 
%
%       - the non-linear response of groundwater head to climate forcing
%       can be accounted for by transforming the input forcing data.
%       Currently, a simple 1-D vertically integrated soil model is
%       available (see climateTransform_soilMoistureModels), allowing the
%       construction of a soil model with between 1 and 7 parameters.
%
%       - Pumping drawdown can be simulated using the Ferris-Knowles 
%       response function for confined aquifers, the Ferris-
%       Knowles-Jacobs function for an unconfined aquifer
%       or the Hantush function for a leaky aquifer. Each of these
%       functions allow multiple production bores to be accounted for and
%       each can have recharge or no-flow boundary conditions.
%
%       - The influence of streamflow can be approximated using the Bruggeman 
%       reponse function. Importantly, this function is a prototype and
%       should be used with caution. It is likely that the covariance
%       between rainfall and streamflow complicates the estimation of the
%       impacts of streamflow.
%
%       - a model can be fit to irregularly sampled groundwater head
%       observations via use of an exponential noise function.
%
%       - new data transfer functions can easily be defined simply
%       by the creation of new response function class definitions.
%
%       - calibration of the model is not to the observed head, but to the
%       innovations between time steps. That is, the residual between the
%       observed and modelled head is derived and the innovation is
%       calculated as the prior residual minus the later residual
%       multiplied by the exponental noise function.
%
%       - the contribution of a given driver to the head is calculated as
%       the integral of the historic forcing times the weighting function.
%       This is undertaken by get_h_star() and the mex .c compiled function
%       doIRFconvolution(). If the forcing is instantaneous (for example a
%       flux a soil moisture model) then Simpon's integration is
%       undertaken. However, if the forcing is a daily integral (such as
%       precipitation or daily pumping volumes) then daily trapazoidal
%       integration of the weighting function is undertaken and then 
%       multiplied by the daily flux.
%
%   Below are links to the details of the public methods of this class. See
%   the 'model_TFN' constructor for details of how to build a model.
%
% See also:
%   HydroSight: time_series_model_calibration_and_construction;
%   model_TFN: model_construction;
%   calibration_finalise: initialisation_of_model_prior_to_calibration;
%   calibration_initialise: initialisation_of_model_prior_to_calibration;
%   get_h_star: main_method_for_calculating_the_head_contributions.
%   getParameters: returns_a_vector_of_parameter_values_and_names;
%   objectiveFunction: returns_a_vector_of_innovation_errors_for_calibration;
%   setParameters: sets_model_parameters_from_input_vector_of_parameter_values;
%   solve: solve_the_model_at_user_input_sime_points;
%
% Dependencies:
%   responseFunction_Pearsons.m
%   responseFunction_Hantush.m
%   responseFunction_Bruggeman.m
%   responseFunction_FerrisKnowles.m
%   responseFunction_Hantush.m
%   responseFunction_FerrisKnowlesJacobs.m
%   derivedweighting_PearsonsNegativeRescaled.m
%   derivedweighting_PearsonsPositiveRescaled.m
%   climateTransform_soilMoistureModels.m
%
% References:
%   von Asmuth J. R., Bierkens M. F. P., Mass K., 2002, Transfer
%   dunction-noise modeling in continuous time using predefined impulse
%   response functions.
%
%   von Asmuth J. R., Mass K., Bakker M., Peterson J., 2008, Modeling time
%   series of ground water head fluctuations subject to multiple stresses.
%   Groundwater, 46(1), pp30-40.
%
%   Peterson and Western (2014), Nonlinear time-series modeling of unconfined
%   groundwater head, Water Resources Research, DOI: 10.1002/2013WR014800
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure Engineering, 
%   The University of Melbourne.
%
% Date:
%   26 Sept 2014
       
    properties
    end

%%  STATIC METHODS        
% Static methods used by the Graphical User Interface to inform the
% user of the available model types. Any new models must be listed here
% in order to be accessable within the GUI.
    methods(Static)
        % Provides a simple description of the function for the user
        % interface.
        function str = description()
           
            str = {['"model_TFN" is a highly flexible nonlinear transfer function noise time-series model. ', ...
                   'It allows the statistical modelling of irregular groundwater head observations by the weighting ', ...
                   'of observed forcing data (e.g. daily rainfall, pumping or landuse change).', ...
                   'The forcing data can also be transformed to account for non-linear processes such as runoff or free drainage.'], ...
                   '', ...
                   'The model has the following additional features:', ...
                   '   - time-series extrapolation and interpolation.', ...
                   '   - decomposition of the hydrograph to individual drivers ', ...
                   '     (e.g. seperation of climate and pumping impacts).', ...
                   '   - decomposition of the hydrograph over time ', ...
                   '     (i.e. decomposition to influence from 1, 2, 5 and 10 years ago).', ...
                   '   - estimation of hydraulic properties (if pumping is simulated).', ...
                   '', ...
                   'For further details of the algorithms use the help menu or see:', ...
                   '', ...
                   '   - Peterson, T. J., and A. W. Western (2014), Nonlinear time-series modeling of unconfined groundwater head, Water Resour. Res., 50, 8330–8355, doi:10.1002/2013WR014800.', ...
                   '', ...
                   '   - Shapoori V., Peterson T. J., Western A. W. and Costelloe J. F., (2015), Top-down groundwater hydrograph time series modeling for climate-pumping decomposition. Hydrogeology Journal'};                     
               
        end        
    end        
    
%%  PUBLIC METHODS              
    methods
        
%% Model constructor
        function obj = model_TFN_HMM(bore_ID, obsHead, forcingData_data,  forcingData_colnames, siteCoordinates, varargin)           
            obj = obj@model_TFN(bore_ID, obsHead, forcingData_data,  forcingData_colnames, siteCoordinates, varargin{1})    %The @ is here for inheritancy, so it will run model_TFN, and then run this

            obj.parameters = rmfield(obj.parameters,'noise');
            obj.parameters.noise1 = noise(obsHead(:,1));
            obj.parameters.noise2 = noise(obsHead(:,1));
            %obj.parameters.noise = rmfield(obj.parameters.noise,'alpha');
            %obj.parameters.noise.alpha_1 = -1;
            %obj.parameters.noise.alpha_2 = -1;
            obj.parameters.datum.d1 = mean(obsHead(:,end));
            obj.parameters.datum.d2 = mean(obsHead(:,end));
            obj.parameters.probabilities.initial_1 = 0.5;
            obj.parameters.probabilities.initial_2 = 0.5;
            obj.parameters.probabilities.trans_1 = 0.5;
            obj.parameters.probabilities.trans_2 = 0.5;

        end
%% Calculate objective function vector. 
        function [objFn_1,objFn_2, h_star1,h_star2, colnames, drainage_elevation] = objectiveFunction(params, time_points, obj, varargin)
% objectiveFunction calculates the objective function vector. 
%
% Syntax:
%   [objFn, h_star, colnames, drainage_elevation] = objectiveFunction(params,time_points, obj)
%
% Description:
%   Solves the model for the input parameters and calculates the objective
%   function vector. Importantly, the objective function vector is not
%   simply the difference between the observed and modelled heads. Because
%   the model uses a noise model, the residual between the observed 
%   and modelled head is first derived and then the innovation
%   is calculated as the prior residual minus the later residual multiplied
%   by the exponental noise function. Finally, the objective function
%   weights this vector according to the time-step between observations.
%
%   Imporantly, the numerator of the weighting equation from von Asmuth et al 2002
%   rounds to zero when the number of samples is very large. This occurs
%   because it is effecively a geometric mean and its product term for n
%   (where n is the number of observation for calibration minus 1) rounds
%   to zero as a result of machine precision. This was overcome by adoption
%   of a restructuring of the geometric meazn in term of exp and log terms. 
%
% Inputs:
%   params - column vector of the optima parameters derived from
%   calibration.
%
%   time_points - column vector of the time points to be simulated.  
%
%   obj -  model object
%
% Outputs:
%   objFn - scalar objective function value.
%
%   h_star - matrix of the contribution from various model components and
%   their summed influence. The matrix columns are in the order of:
%   date/time, summed contribution to the head, contribution from
%   component i.
%
%   colnames - column names for matrix 'head'.
%
%   drainage_elevation - drainage elevation constant.
%
% Example:
%   see HydroSight: time_series_model_calibration_and_construction;
%
% See also:
%   HydroSight: time_series_model_calibration_and_construction;
%   model_TFN: model_construction;
%   calibration_finalise: initialisation_of_model_prior_to_calibration;
%   calibration_initialise: initialisation_of_model_prior_to_calibration;
%   get_h_star: main_method_for_calculating_the_head_contributions.
%   getParameters: returns_a_vector_of_parameter_values_and_names;
%   setParameters: sets_model_parameters_from_input_vector_of_parameter_values;
%   solve: solve_the_model_at_user_input_sime_points;
%
% References:
%   von Asmuth J. R., Bierkens M. F. P., Mass K., 2002, Transfer
%   dunction-noise modeling in continuous time using predefined impulse
%   response functions.
%
% Author: 
%   Dr. Tim Peterson, The Department of Infrastructure
%   Engineering, The University of Melbourne.
%
% Date:
%   26 Sept 2014    
            
            % Set model parameters
            setParameters(obj, params, obj.variables.param_names);
            
            % If varargin is a structural variable then the model
            % simulation is to use predefined values of the drainage
            % elevation and the mean forcing. Note, these inputs are only
            % to be provided if not doing simulation.
            getLikelihood = false;
            drainage_elevation=[];
            mean_forcing=[];
            if ~isempty(varargin)
                if isstruct(varargin{1})
                    drainage_elevation=varargin{1}.drainage_elevation;
                    mean_forcing=varargin{1}.mean_forcing;
                elseif islogical(varargin{1})
                    getLikelihood=varargin{1};
                else
                    error('The input varargin must be either a logical or a structural variable.');
                end
            end
            
            % Calc deterministic component of the head.
            if isempty(mean_forcing)
                [h_star, colnames] = get_h_star(obj, time_points);                 
            else
                [h_star, colnames] = get_h_star(obj, time_points,mean_forcing);                 
            end

            % Increment count of function calls
            obj.variables.nobjectiveFunction_calls = obj.variables.nobjectiveFunction_calls +  size(params,2);
                        
            % Return of there are nan or inf value
            if any(isnan(h_star(:,2)) | isinf(h_star(:,2)))
                if getLikelihood
                    objFn_1 = -inf;
                    objFn_2 = -inf;
                else
                    objFn_1 = inf;
                    objFn_2 = inf;
                end
                return;
            end
               

            h_star1 = h_star;
            h_star1(:,2) = h_star(:,2) + obj.parameters.datum.d1;
            h_star2 = h_star;
            h_star2(:,2) = h_star(:,2) + obj.parameters.datum.d2;
            
            % If the results from this method call are not to be used for
            % summarising calibration results, then exit here. This is
            % required because the innovations can only be calculated at
            % time points for which there are observations. 
            if ~obj.variables.doingCalibration
                objFn_1 = [];
                objFn_2 = [];
                return;
            end
            
            % Calculate residual between observed and modelled. 
            t_filt = find( obj.inputData.head(:,1) >=time_points(1)  ...
                & obj.inputData.head(:,1) <= time_points(end) );          
            resid1= obj.inputData.head(t_filt,2)  - h_star1(:,2);
            resid2= obj.inputData.head(t_filt,2)  - h_star2(:,2);
            %REPLACE ALPHA_1 with OBJECT.NOISE
            
            % Calculate innovations using residuals from the deterministic components.            
            innov1 = resid1(2:end) - resid1(1:end-1).*exp( -10.^obj.parameters.noise.alpha_1 .* obj.variables.delta_time );
            innov2 = resid2(2:end) - resid2(1:end-1).*exp( -10.^obj.parameters.noise.alpha_2 .* obj.variables.delta_time );
            
            % Calculate objective function
            objFn_1 = exp(mean(log( 1- exp( -2.*10.^obj.parameters.noise.alpha_1 .* obj.variables.delta_time) ))) ...
                    ./(1- exp( -2.*10.^obj.parameters.noise.alpha_1 .* obj.variables.delta_time )) .* innov1.^2;
            objFn_2 = exp(mean(log( 1- exp( -2.*10.^obj.parameters.noise.alpha_2 .* obj.variables.delta_time) ))) ...
                    ./(1- exp( -2.*10.^obj.parameters.noise.alpha_2 .* obj.variables.delta_time )) .* innov2.^2;    %Von Asmuth 2005 Paper, see Tim's 2014 paper
  
%%CYCLE THROUGH FOR LOOP EACH TIMESTEP FOR THE LENGTH OF objFN
            % Calculate log liklihood    
            if getLikelihood
                N1 = size(resid1,1);
                N2 = size(resid2,1);
                objFn_1 = -0.5 * N1 * ( log(2*pi) + log(objFn_1./N1)+1);
                objFn_2 = -0.5 * N2 * ( log(2*pi) + log(objFn_2./N2)+1);
            end
        end
        function [timeseries,integers] = viterbi(obj)
            timeseries = ;
            integers = ;    %1 = state 1, 2 = state 2
        end
    end
end