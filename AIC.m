close all; clear; clc;
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/cleaned/");

filename = "47996";

load(filename + "_model_6params_1.mat");
loglikelihood = model_6params.model.variables.objFn;
parameters_number = length(model_6params.model.variables.param_names);
aic = 2 * parameters_number + 2 * loglikelihood;

load(filename + "_model_6params_2.mat");
loglikelihood_2 = model_6params.model.variables.objFn;
parameters_number_2 = length(model_6params.model.variables.param_names);
aic_2 = 2 * parameters_number_2 + 2 * loglikelihood_2;