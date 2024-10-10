classdef model_TFN_HMM_2 < model_TFN
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
        function obj = model_TFN_HMM_2(bore_ID, obsHead, forcingData_data,  forcingData_colnames, siteCoordinates, varargin)
            obj = obj@model_TFN(bore_ID, obsHead, forcingData_data,  forcingData_colnames, siteCoordinates, varargin{1})    %The @ is here for inheritancy, so it will run model_TFN, and then run this

            obj.parameters = rmfield(obj.parameters,'noise');
            obj.parameters.noise1 = noise(obsHead(:,1));
            obj.parameters.datum1 = datum(obsHead);

        end
        %% Calculate objective function vector.
        function [objFn, h_star1, colnames] = objectiveFunction(params, time_points, obj, varargin)
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

            % % If varargin is a structural variable then the model
            % % simulation is to use predefined values of the drainage
            % % elevation and the mean forcing. Note, these inputs are only
            % % to be provided if not doing simulation.
            % getLikelihood = false;
            %drainage_elevation=[];
            mean_forcing=[];
            getLikelihood = false;
            if ~isempty(varargin)
                if isstruct(varargin{1})
                    %drainage_elevation=varargin{1}.drainage_elevation;
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
                [h_star, colnames] = get_h_star(obj, time_points, mean_forcing);
            end

            % if isempty(drainage_elevation)
            %     drainage_elevation = obj.variables.h_bar - mean(h_star(:,2));
            % end

            % Increment count of function calls
            obj.variables.nobjectiveFunction_calls = obj.variables.nobjectiveFunction_calls +  size(params,2);

            % Return of there are nan or inf value
            if any(isnan(h_star(:,2)) | isinf(h_star(:,2)))
                if getLikelihood
                    objFn = -inf;
                else
                    objFn = inf;
                end
                return;
            end

            h_star1 = h_star;
            h_star1(:,2) = h_star(:,2) + obj.parameters.datum1.d;


            
            % if isempty(drainage_elevation)
            %     drainage_elevation = obj.variables.h_bar - mean(h_star(:,2));
            % end

            % If the results from this method call are not to be used for
            % summarising calibration results, then exit here. This is
            % required because the innovations can only be calculated at
            % time points for which there are observations.
            if ~obj.variables.doingCalibration
                objFn = [];
                %objFn_1 = [];
                %objFn_2 = [];
                return;
            end

            %Calculate residual between observed and modelled.
            t_filt = find( obj.inputData.head(:,1) >=time_points(1)  ...
                & obj.inputData.head(:,1) <= time_points(end) );
            resid1= obj.inputData.head(t_filt,2)  - h_star1(:,2);
            delta_time = [inf;obj.variables.delta_time];
            % Calculate innovations using residuals from the deterministic components.
            innov1 = resid1(2:end) - resid1(1:end-1).*exp( -10.^obj.parameters.noise1.alpha .* delta_time(2:end) );
            % at first timestep prior contributions is zero
            innov1 = [resid1(1);innov1];
            

            obj.parameters.noise1.sigma_n = sqrt(mean( innov1.^2 ./ (1 - exp( -2 .* 10.^obj.parameters.noise1.alpha .* delta_time ))));
            objFn = -sum(log(pdf('Normal', resid1, 0, obj.parameters.noise1.sigma_n)));
                
        end

        %% Finalise the model following calibration.
        function h_star = calibration_finalise(obj, params, useLikelihood)
            % calibration_finalise finalises the model following calibration.
            %
            % Syntax:
            %   calibration_finalise(obj, params)
            %
            % Description:
            %   Finalises the model following calibration and assigns the final
            %   parameters and additional variables to the object for later simulation.
            %   Of the variables calculated, the most essential for the model
            %   is the scalar drainage, obj.variables.d. Other variables that are also
            %   important include:
            %       - a vector of innovations,  obj.variables.innov, for detection
            %         of serial correlation in the model errors;
            %       - the noise standard deviation, obj.variables.sigma_n.
            %
            % Input:
            %   obj -  model object
            %
            %   params - column vector of the optima parameters derived from
            %   calibration.
            %
            % Outputs:
            %   h_star: simulated head is returned.
            %
            % Example:
            %   see HydroSight: time_series_model_calibration_and_construction;
            %
            % See also:
            %   HydroSight: time_series_model_calibration_and_construction;
            %   model_TFN: model_construction;
            %   calibration_initialise: initialisation_of_model_prior_to_calibration;
            %   get_h_star: main_method_for_calculating_the_head_contributions.
            %   getParameters: returns_a_vector_of_parameter_values_and_names;
            %   objectiveFunction: returns_a_vector_of_innovation_errors_for_calibration;
            %   setParameters: sets_model_parameters_from_input_vector_of_parameter_values;
            %   solve: solve_the_model_at_user_input_sime_points;
            %
            % Author:
            %   Dr. Tim Peterson, The Department of Infrastructure
            %   Engineering, The University of Melbourne.
            %
            % Date:
            %   26 Sept 2014

            % Set the stochastic forcing to NOT be in a 'calibration' state.
            setStochForcingState(obj, false, obj.variables.t_start, obj.variables.t_end);

            % Initialise data
            nparamSets = size(params,2);
            objFn = NaN(1,nparamSets);
            time_points = obj.variables.time_points;
            companantNames = fieldnames(obj.inputData.componentData);
            nCompanants = size(companantNames,1);
            forcingMean = inf(nCompanants,1,nparamSets);
            nvars = zeros(nCompanants,1);

            % Run objective function to get data and num cols of h_star =
            % in case there's >1 parameter set.
            [objFn(:,1), h_star] = objectiveFunction(params(:,1), time_points, obj,useLikelihood);

            t_filt = find( obj.inputData.head(:,1) >=time_points(1)  ...
                & obj.inputData.head(:,1) <= time_points(end) );

            for i=1:nparamSets
            %Calculate residual between observed and modelled.
                resid1= obj.inputData.head(t_filt,2)  - h_star(:,2,i);

                delta_time = [inf;obj.variables.delta_time];

                innov1 = resid1(2:end) - resid1(1:end-1).*exp( -10.^obj.parameters.noise1.alpha(i) .* delta_time(2:end) );
                innov1 = [resid1(1);innov1];
                
                obj.parameters.noise1.sigma_n(i) = sqrt(mean( innov1.^2 ./ (1 - exp( -2 .* 10.^obj.parameters.noise1.alpha(i) .* delta_time ))));
            end

            noise = getNoise(obj.parameters.noise1, time_points);

            h_star = [h_star(:,:,:), h_star(:,2,:) - noise(:,2,:), ...
                h_star(:,2,:) + noise(:,3,:)];

            for j=1:nCompanants
                nvars(j) = size(obj.variables.(companantNames{j}).forcingData,2);
                forcingMean(j,1:nvars(j),1) = mean(obj.variables.(companantNames{j}).forcingData,1);
            end

            for j=1:nCompanants
                obj.variables.(companantNames{j}).forcingMean = forcingMean(j,1:nvars(j),:);
            end

            obj.variables.residualMean = mean(abs(resid1));
            obj.variables.d = obj.parameters.datum1.d;
            obj.variables.objFn = objFn;
            clear d objFn forcingMean


            % Set a flag to indicate that calibration is complete.
            obj.variables.doingCalibration = false;

            % Free memory within mex function
            try
                junk=doIRFconvolutionPhi([], [], [], [], false, 0);            %#ok<NASGU>
            catch
                % continue
            end
        end
        %% Solve Function
        function [h_star, colnames, noise] = solve(obj, time_points)
            [params, param_names] = getParameters(obj);
            nparamsets = size(params,2);

            companants = fieldnames(obj.inputData.componentData);
            nCompanants = size(companants,1); 
            for ii=1:nparamsets
                % Get the calibration estimate of the mean forcing for the
                % current parameter set. This is a bit of a work around to
                % handle the issue of each parameter set having a unique
                % mean forcing (if a forcing transform is undertaken). The
                % workaround was required when DREAM was addded.
                for j=1:nCompanants                    
                    calibData(ii,1).mean_forcing.(companants{j}) = obj.variables.(companants{j}).forcingMean(:,:,ii); %#ok<AGROW> 
                end                
            end
            %ERROR, HEAD IS COMING OUT AS A SINGLE VALUE REPEATED
            [~, h_star,colnames] = objectiveFunction(params(:,1), time_points, obj, calibData(1)); %Check ~
            noise = getNoise(obj.parameters.noise1, time_points);
           
            
            figure()
            head = getObservedHead(obj);
            t = head(:,1);
            plot(t, head(:,2),'.-k', 'LineWidth',0.025,'MarkerSize',2 );
           
            hold on
            plot(t, h_star(:,2),'.-b','LineWidth',0.025); 
            legend("Observed","Simulated");

            tspan = min(t):1:max(t);
            filt = day(tspan)==1 & month(tspan)==1;
            tspan = tspan(filt);
            yspan = year(tspan);
            xticks(gca,tspan);
            xticklabels(gca,yspan);
            hold off

            if nparamsets>1
                h_star = cat(3, h_star, zeros(size(h_star,1),size(h_star,2), nparamsets-1));
                parfor jj=2:size(params,2)
                    [~, h_star(:,:,jj)] = objectiveFunction(params(:,jj), time_points, obj, calibData(jj));
                end
            end

            % Set the parameters if >1 parameter sets
            if nparamsets>1
                setParameters(obj,params, param_names);
            end

            % Clear matrix of indexes for tor at each time_points
            obj.variables = rmfield(obj.variables, 'theta_est_indexes_min');
            obj.variables = rmfield(obj.variables, 'theta_est_indexes_max');

        end
    end
end



