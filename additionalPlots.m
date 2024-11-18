close all; clear; clc;
addpath(genpath("downloads"));
addpath("downloads/models/");
set(groot,'defaultfigureposition',[0,0,2560,1440]);
list = {'46854','47996','54925','56252','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';
fourPlotsOneAxes = false;
twoStatePlot = false;
oneStatePlot = false;
vals = 2; %1:length(list);

if fourPlotsOneAxes == true
    for k = vals %1:length(list)
        filename = list{k};

        % Load Single State model
        load(filename + "_model_6params_1.mat");
        time_obs = model_6params.model.inputData.head(:,1);
        head_obs = model_6params.model.inputData.head(:,2);
        time_steps_single = model_6params.simulationResults{1,1}.head(:,1);
        head_1 = model_6params.simulationResults{1,1}.head(:,2);

        clear model_6params;

        % Load Two State model
        load(filename + "_model_6params_2.mat");
        time_steps_two = model_6params.simulationResults{1,1}.head(:,1);
        head_2 = model_6params.simulationResults{1,1}.head(:,2);

        datum1 = model_6params.model.parameters.datum1;
        datum2 = model_6params.model.parameters.datum2;

        if datum1.d > datum2.d
            wetState = 1;
            dryState = 2;
        elseif datum2.d > datum1.d
            wetState = 2;
            dryState = 1;
        else
            error("Datum 1 or 2 either aren't defined or are equal");
        end

        figure(k);
        plot(time_steps_single,head_1,'color', [.5 .5 .5],'LineStyle','--','LineWidth',2);
        hold on;
        plot(time_obs,head_obs,'k','LineWidth',2);

        stateHead = model_6params.model.variables.StateHead;
        stateHead_1 = stateHead;
        stateHead_2 = stateHead;

        for j = 1:size(stateHead_1, 1)
            if stateHead_2(j, 3) == wetState
                stateHead_2(j, 2) = NaN;
            end
        end

        t = time_obs;
        plot(t, stateHead(:,2),'.-b','LineWidth',2);
        % plot(t,stateHead_1(:,2),"Marker",".","Color","b");
        plot(t,stateHead_2(:,2),"Marker",".","Color",'r','LineWidth',2);
        tspan = min(t):1:max(t);
        filt = day(tspan)==1 & month(tspan)==1;
        tspan = tspan(filt);
        yspan = year(tspan);
        xticks(gca,tspan);
        xticklabels(gca,yspan);
        title(['ID' filename ' Two-State Solved Plot Compared to Single-State Solved Plot']);
        legend('Single State Model','Observed Groundwater Head','Wet State','Dry State');
        ylabel('Groundwater Head (m)');
        fontsize(36,"points");
        hold off
        saveas(gca,['downloads/plots_final/' filename '_solvedPlot.png']);
        close all;

        % ADD Additional h_stars from model_6params.model.variables.h_stars
    end
end

if fourPlotsOneAxes == false && twoStatePlot == false && oneStatePlot == false
    for k = vals
        filename = list{k};

        % Load Two State model
        load(filename + "_model_6params_2.mat");
        time_steps = model_6params.simulationResults{1,1}.head(:,1);
        head_obs = model_6params.model.inputData.head(:,2);
        head_model1 = model_6params.model.variables.h_stars(:,2); %hstars(:,2)
        head_model2 = model_6params.model.variables.h_stars(:,3); %hstars(:,3)

        figure(k);
        plot(time_steps,head_obs,'k','LineWidth',2);
        hold on;

        dryIndex = model_6params.model.variables.StateHead;
        plot(time_steps, dryIndex(:,2),'.-b','LineWidth',2);

        for j = 1:size(dryIndex, 1)
            if dryIndex(j, 3) == 1
                dryIndex(j, 2) = NaN;
            end
        end
        plot(time_steps(~isnan(dryIndex(:,2))), head_model1(~isnan(dryIndex(:,2))),'.-','Color',[0.3010 0.7450 0.9330],'LineWidth',2);

        plot(time_steps,dryIndex(:,2),"Marker",".","Color","r",'LineWidth',2);
        tspan = min(time_steps):1:max(time_steps);
        filt = day(tspan)==1 & month(tspan)==1;
        tspan = tspan(filt);
        yspan = year(tspan);
        xticks(gca,tspan);
        xticklabels(gca,yspan);
        title(['ID' filename ' Two-State Solved Plot Comparing States After Switch']);
        legend('Observed Groundwater Head','Wet State','Wet State after Switch','Dry State');
        fontsize(20,"points");
        ylabel('Groundwater Head (m)');
        hold off

        % ADD Additional h_stars from model_6params.model.variables.h_stars
    end
end

if twoStatePlot == true
    for k = vals %1:length(list)
        filename = list{k};

        % Load Single State model
        load(filename + "_model_6params_1.mat");
        time_obs = model_6params.model.inputData.head(:,1);
        head_obs = model_6params.model.inputData.head(:,2);
        time_steps_single = model_6params.simulationResults{1,1}.head(:,1);

        clear model_6params;

        % Load Two State model
        load(filename + "_model_6params_2.mat");
        time_steps_two = model_6params.simulationResults{1,1}.head(:,1);
        head_2 = model_6params.simulationResults{1,1}.head(:,2);

        datum1 = model_6params.model.parameters.datum1;
        datum2 = model_6params.model.parameters.datum2;

        if datum1.d > datum2.d
            wetState = 1;
            dryState = 2;
        elseif datum2.d > datum1.d
            wetState = 2;
            dryState = 1;
        else
            error("Datum 1 or 2 either aren't defined or are equal");
        end

        figure(k);
        plot(time_obs,head_obs,'k','LineWidth',2);
        hold on;

        stateHead = model_6params.model.variables.StateHead;
        stateHead_1 = stateHead;
        stateHead_2 = stateHead;

        for j = 1:size(stateHead_1, 1)
            if stateHead_2(j, 3) == wetState
                stateHead_2(j, 2) = NaN;
            end
        end

        t = time_obs;
        plot(t, stateHead(:,2),'.-b','LineWidth',2);
        plot(t,stateHead_2(:,2),"Marker",".","Color",'r','LineWidth',2);
        tspan = min(t):1:max(t);
        filt = day(tspan)==1 & month(tspan)==1;
        tspan = tspan(filt);
        yspan = year(tspan);
        xticks(gca,tspan);
        xticklabels(gca,yspan);
        title(['ID' filename ' Two-State Solved Plot']);
        legend('Observed Groundwater Head','Wet State','Dry State');
        fontsize(36,"points");
        hold off

        % ADD Additional h_stars from model_6params.model.variables.h_stars
    end
end

if twoStatePlot == false && oneStatePlot == true
    for k = vals %1:length(list)
        filename = list{k};

        % Load Single State model
        load(filename + "_model_6params_1.mat");
        time_obs = model_6params.model.inputData.head(:,1);
        head_obs = model_6params.model.inputData.head(:,2);
        time_steps_single = model_6params.simulationResults{1,1}.head(:,1);
        head_1 = model_6params.simulationResults{1,1}.head(:,2);

        clear model_6params;

        figure(k);
        plot(time_steps_single,head_1,'b','LineWidth',2);
        hold on;
        plot(time_obs,head_obs,'k','LineWidth',2);
        t = time_obs;
        tspan = min(t):1:max(t);
        filt = day(tspan)==1 & month(tspan)==1;
        tspan = tspan(filt);
        yspan = year(tspan);
        xticks(gca,tspan);
        xticklabels(gca,yspan);
        title(['ID' filename ' One-State Solved Plot']);
        legend('Single State Model','Observed Groundwater Head');
        fontsize(36,"points");
        hold off

        % ADD Additional h_stars from model_6params.model.variables.h_stars
    end
end