close all; clear; clc;
set(groot,'defaultfigureposition',[0,0,2560,1440]);
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/cleaned/");
addpath("downloads/plots/");

list = {'46854','47996','54925','56252','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';
% 47996, 82095, 105287, 108944, 110197ii, 111691, 112236, 119429
% i = [2, 12, 17, 18, 21, 25, 26, 27];

load("forcing_northingeasting.mat");
load("bore_northingeasting.mat");

residualMean = zeros(length(list),2);

two_stateMethod = true;
single_stateMethod = true;

for i = 2 %1:length(list)

    filename = list{i};

    load(filename + "_boreData.mat");
    load(filename + "_forcingData.mat");

    % Reformat the matric of forcing data to a sturctire variable containing
    % the column names.
    forcingDataStruct.data = forcingData;
    forcingDataStruct.colnames = {'YEAR','MONTH','DAY','PRECIP','APET','RevegFrac'};

    bore_ID = convertStringsToChars("ID" + filename);
    bore_n = bore_ne(i,2);
    bore_e = bore_ne(i,3);
    forcing_n = forcing_northingeasting(i,2);
    forcing_e = forcing_northingeasting(i,3);

    siteCoordinates = {bore_ID, bore_n, bore_e; 'PRECIP', forcing_n, forcing_e; 'APET', forcing_n, forcing_e; 'RevegFrac',forcing_n, forcing_e};

    forcingTransform_Precip = {'transformfunction', 'climateTransform_soilMoistureModels'; ...
        'forcingdata', {'precip','PRECIP';'et','APET'}; ...
        'outputdata', 'drainage'; ...
        'options', {'SMSC',2,[];'beta',0,'';'k_sat',1,'fixed';'alpha',0,'fixed'}};

    forcingTransform_ET = {'transformfunction', 'climateTransform_soilMoistureModels_v2'; ...
        'outputdata', 'evap_gw_potential'};

    modelOptions_6params = { 'precip','weightingfunction','responseFunction_Pearsons'; ...
        ['precip' ...
        ''],'forcingdata',forcingTransform_Precip};
    maxObsFreq = 7;
    % Build the 6 parameter model.
    modelLabel = bore_ID;

    if two_stateMethod == true
        % HMM 2 State Model
        model_6params = HydroSightModel(modelLabel, bore_ID, 'model_TFN_HMM', boreDataWL, maxObsFreq, forcingDataStruct, siteCoordinates, modelOptions_6params);
        SchemeSetting.ngs = 12*12;
        calibrateModel(model_6params, [], 0, inf, 'SP-UCI', SchemeSetting);
        calibrateModelPlotResults(model_6params,[]);
        calibratedPlot_twoState = gca;
        saveas(gca,['downloads/plots/' filename '_calibratedPlot.png']);

        time_points = model_6params.model.variables.time_points;
        newForcingData = [];
        simulationLabel = 'default simulation';
        doKrigingOnResiduals = false;

        solveModel(model_6params, time_points, newForcingData, simulationLabel, doKrigingOnResiduals);
        title([bore_ID ' Two-State Solved Plot']);
        solvedPlot_twoState = gca;
        saveas(gca,['downloads/plots/' filename '_solved2StatePlot.png']);
        residualMean(i,1) = model_6params.model.variables.residualMean;
        save(['downloads/plots/' filename '_model_6params_2.mat'],"model_6params");
        clear model_6params;
    end

    if single_stateMethod == true
        % HMM 1 State Model
        model_6params = HydroSightModel(modelLabel, bore_ID, 'model_TFN_HMM_2', boreDataWL, maxObsFreq, forcingDataStruct, siteCoordinates, modelOptions_6params);
        SchemeSetting.ngs = 12*12;
        calibrateModel(model_6params, [], 0, inf, 'SP-UCI', SchemeSetting);
        calibrateModelPlotResults(model_6params,[]);
        calibratedPlot_singleState = gca;
        saveas(gca,['downloads/plots/' filename '_calibratedPlot_singlestate.png']);

        time_points = model_6params.model.variables.time_points;
        newForcingData = [];
        simulationLabel = 'default simulation';
        doKrigingOnResiduals = false;

        solveModel(model_6params, time_points, newForcingData, simulationLabel, doKrigingOnResiduals);
        title([bore_ID ' One-State Solved Plot']);
        solvedPlot_singleState = gca;
        saveas(gca,['downloads/plots/' filename '_solved1StatePlot.png']);
        residualMean(i,2) = model_6params.model.variables.residualMean;
        save(['downloads/plots/' filename '_model_6params_1.mat'],"model_6params");
    end
    close all;
end

% Save a zip file of the plots as a back up that can more easily be moved
% zip('currentPlots.zip','downloads\plots');