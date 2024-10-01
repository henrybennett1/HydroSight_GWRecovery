close all; clear; clc;
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/cleaned/");

list = {'46854','47996','54925','56252','62427','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';

for i = 1:length(list)
    filename = list{i};

    load(filename + "_boreData.mat");
    figure(i);
    plot(boreDataWL(:,5));
end