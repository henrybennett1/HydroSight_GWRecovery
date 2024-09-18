% example_TFN_model:
%
% Description
%   This example builds and calibrates a nonlinear transfer-function noise
%   model. The example is taken from Peterson & Western (2014). The model
%   requires the following three data files: 124705_boreData.mat,
%   124676_boreData.mat and 124705_forcingData.mat.
%
%   By default the example models bore 124676. By commenting out line 29
%   and un-commenting line 28, bore 124705 can be modelled.
%
%   Also, logical variables at line 112 and 113 also define which of 
%   two model structures are to be calibrated.
%
% References:
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
%

%clear all

% Comment out the one bore ID that you DO NOT want to model.
bore_ID = 'ID124705';

if strcmp(bore_ID,'ID124705')
    load('124705_boreData.mat');    %   HydroSight_GWRecovery\algorithms\models\TransferNoise\Example_model\
else
    load('124676_boreData.mat');
end

load('124705_forcingData.mat');



% Reformat the matric of forcing data to a sturctire variable containing
% the column names.
forcingDataStruct.data = forcingData;
forcingDataStruct.colnames = {'YEAR','MONTH','DAY','PRECIP','APET','RevegFrac'};

% To increase performance, we can reduce the length of the climate record.
% This may cause longer time scales to be less reliably estimated.
yearsOfPriorForcing = 100;
forcingData_thresholddate  = datenum( boreDataWL(1,1)- yearsOfPriorForcing, boreDataWL(1,2), boreDataWL(1,3)); 
filt = datenum(forcingDataStruct.data(:,1), forcingDataStruct.data(:,2), forcingDataStruct.data(:,3)) >= forcingData_thresholddate;
forcingDataStruct.data = forcingDataStruct.data(filt,:);

% Define the bore ID and create sume dummy site coordinates. This must be
% for the bore and each column in the forcing file.
siteCoordinates = {bore_ID, 100, 100; 'PRECIP', 100, 100; 'APET', 100, 100; 'RevegFrac',602, 100};

% Define the way in which the precipitation is transformed. In this case it
% is transformed using the 'climateTransform_soilMoistureModels' soil
% model. 
% Next, the soil ODE needs inputs data 'precip' and 'et' and the forcing
% data input columns are 'PRECIP' and 'ET'.
% Next, the 'outputdata' that is to be taken from the soil model is
% defined. Each model has fixed options and here we're taking
% 'drainage'.
% Lastly, we can set 'options' for the soil model. In this case we are
% defining the initial values for three parameters (SMSC, beta, ksat) and
% fixing alpha to zero.
forcingTransform_Precip = {'transformfunction', 'climateTransform_soilMoistureModels'; ...
               'forcingdata', {'precip','PRECIP';'et','APET'}; ...
               'outputdata', 'drainage'; ...
               'options', {'SMSC',2,[];'beta',0,'';'k_sat',1,'fixed';'alpha',0,'fixed'}};
           
% The transformation of the ET is then defined. However because we've already 
% defined the soil model, we only need to specify the output we require.
% Here we're selecting  'evap_gw_potential', which is the potential ET -
% actual soil ET.
forcingTransform_ET = {'transformfunction', 'climateTransform_soilMoistureModels_v2'; ...
               'outputdata', 'evap_gw_potential'};            

modelOptions_6params = { 'precip','weightingfunction','responseFunction_Pearsons'; ...
                        ['precip' ...
                        ''],'forcingdata',forcingTransform_Precip};                    
                    
% Set the maximum frequency of water level obs
maxObsFreq = 7;

% Select which model structures to build and calibrate.
run6paramModel = true;

% Define a model label
modelLabel = 'Great Western Catchment - no landuse change';

if run6paramModel
    % Build the 6 parameter model.
    model_6params = HydroSightModel(modelLabel, bore_ID, 'model_TFN_HMM', boreDataWL, maxObsFreq, forcingDataStruct, siteCoordinates, modelOptions_6params);

    % Set the number of SP-UCI calibration clusters per parameter
    SchemeSetting.ngs = 4*12;
    
    % Calibrate the 6 parameter model.
    calibrateModel(model_6params, [], 0, inf, 'SP-UCI', SchemeSetting);
    %calibrateModel(model_6params, [], 0, inf, 'CMA-ES', SchemeSetting);
    
    % Plot the calibration results.    
    calibrateModelPlotResults(model_6params,[]);
   
    % Plot the simulation results. 
    time_points = model_6params.model.variables.time_points;
    newForcingData = [];
    simulationLabel = 'default simulation';
    doKrigingOnResiduals = false;    
    
    solveModel(model_6params, time_points, newForcingData, simulationLabel, doKrigingOnResiduals);    
    solveModelPlotResults(model_6params, simulationLabel, []);  
end

figure(2);
hold on;
plot(model_6params.model.variables.time_points,model_6params.model.variables.viterbiStates/4+267.15,'r')