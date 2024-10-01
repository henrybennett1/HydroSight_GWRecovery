close all; clear; clc;
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/cleaned/");

list = {'46854','47996','54925','56252','62427','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';

coords = NaN(length(list),3);
for i = 1:length(list)

    filename = list{i};
    boreDataWL = readtable(['downloads/raw/' filename '_boreData.xlsx'],'ReadVariableNames',false); 
    forcingData = readtable(['downloads/raw/' filename '_ForcingData.csv'],'ReadVariableNames',false);
    % Remove non AHD values
    boreFilt = strcmp(boreDataWL{:,"Var6"}, 'Bore Water Level AHD');
    boreDataWL = boreDataWL(boreFilt,:);

    % Remove 0 values
    boreFilt = logical(boreDataWL{:,"Var8"});
    boreDataWL = boreDataWL(boreFilt,:);

    % Remove negative values
    boreFilt = logical((sign(boreDataWL{:,"Var8"})+1)/2);
    boreDataWL = boreDataWL(boreFilt,:);

    % Remove 99th percentile outliers
    sigma = std(boreDataWL{:,"Var8"});
    mu = mean(boreDataWL{:,"Var8"});
    perc_99 = norminv(0.99,mu,sigma);
    perc_01 = norminv(0.01,mu,sigma);
   
    boreFilt = boreDataWL{:,"Var8"} < perc_99;
    boreDataWL = boreDataWL(boreFilt,:);
    boreFilt = boreDataWL{:,"Var8"} > perc_01;
    boreDataWL = boreDataWL(boreFilt,:);

    % Defining the coordinates
    coords(i,:) = [str2double(filename), forcingData{1,1}, forcingData{1,2}];

    % This could probably be more 'elegant', but I'm just making it work

    temp_bore_datetime = datevec(table2array(boreDataWL(:,{'Var3'})));
    boreDataWL = [temp_bore_datetime(:,1:3), zeros(length(temp_bore_datetime(:,1)),1), ...
        table2array(boreDataWL(:,{'Var8'}))];

    temp_force_datetime = datevec(table2array(forcingData(:,{'Var3'})));
    forcingData = [temp_force_datetime(:,1:3), table2array(forcingData(:,{'Var4'})), ...
        table2array(forcingData(:,{'Var6'})), zeros(length(temp_force_datetime(:,1)),1)];    

    save(['downloads/cleaned/' filename '_boreData.mat'],"boreDataWL");
    save(['downloads/cleaned/' filename '_forcingData.mat'],"forcingData");

    
    
end

save('downloads/cleaned/forcing_coords.mat',"coords");

% forcing data is in format [yy, mm, dd, precip, apet, revegfrac]
% bore data is in format [yy, mm, dd, something, value]