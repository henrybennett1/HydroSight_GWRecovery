close all; clear; clc;
set(groot,'defaultfigureposition',[0,0,1280,720]);
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/cleaned/");
addpath("downloads/plots/");

list = {'46854','47996','54925','56252','62427','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';

load("forcing_coords.mat");

residualMean = zeros(length(list),2);

for i = 1:length(list)
    filename = list{i};

    load(filename + "_boreData.mat");
    load(filename + "_forcingData.mat");

    % Reformat the matric of forcing data to a sturctire variable containing
    % the column names.
    forcingDataStruct.data = forcingData;
    forcingDataStruct.colnames = {'YEAR','MONTH','DAY','PRECIP','APET','RevegFrac'};

    bore_ID = convertStringsToChars("ID" + filename);

    forcing_longitude = coords(i,2);
    forcing_latitude = coords(i,3);
    
    siteCoordinates = {bore_ID, 100, 100; 'PRECIP', 100, 100; 'APET', 100, 100; 'RevegFrac',602, 100}; %FIX
    
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
    
    % HMM 2 State Model
    model_6params = HydroSightModel(modelLabel, bore_ID, 'model_TFN_HMM', boreDataWL, maxObsFreq, forcingDataStruct, siteCoordinates, modelOptions_6params);
    SchemeSetting.ngs = 4*12;
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
    clear model_6params;

    % HMM 1 State Model
    model_6params = HydroSightModel(modelLabel, bore_ID, 'model_TFN_HMM_2', boreDataWL, maxObsFreq, forcingDataStruct, siteCoordinates, modelOptions_6params);
    SchemeSetting.ngs = 4*12;
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
end

save('downloads/plots/meanresiudals.mat',"residualMean");
% Save a zip file of the plots as a back up that can more easily be moved
zip('downloads\currentPlots.zip','downloads\plots');