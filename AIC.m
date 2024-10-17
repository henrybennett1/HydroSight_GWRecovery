close all; clear; clc;
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/models/");

%list = {'46854','47996','62427', '70921', '95076', '98865', '104930', '110738', '111525', '111530', '141235'};
list = {'46854','47996','54925','56252','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';

AIC_values = NaN(length(list),3);

for i = 1:length(list)

    filename = convertCharsToStrings(list(i));

    load(filename + "_model_6params_1.mat");
    loglikelihood = model_6params.model.variables.objFn;
    parameters_number = length(model_6params.model.variables.param_names);
    aic = 2 * parameters_number + 2 * loglikelihood;

    load(filename + "_model_6params_2.mat");
    loglikelihood_2 = model_6params.model.variables.objFn;
    parameters_number_2 = length(model_6params.model.variables.param_names);
    aic_2 = 2 * parameters_number_2 + 2 * loglikelihood_2;
    AIC_values(i,1) = filename;
    AIC_values(i,2) = aic;
    AIC_values(i,3) = aic_2;
end