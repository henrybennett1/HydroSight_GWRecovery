close all; clear; clc;
addpath(genpath("downloads"));
addpath("downloads/raw/");
addpath("downloads/cleaned/");

list = {'46854','47996','54925','56252','62427','63740','64139','66622', ...
    '70921','75563','81957','82095','95076','98865','103344','104930','105287', ...
    '108944','110104','110166','110197','110738','111525','111530','111691','112236', ...
    '119429','141235'}';

fourPlotsOneAxes = true;

if fourPlotsOneAxes == true
    for k = 5 %1:length(list)
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

        figure(k);
        plot(time_steps_single,head_1,'color', [.5 .5 .5],'LineStyle','--');
        hold on;
        plot(time_obs,head_obs,'k');

        stateHead = model_6params.model.variables.StateHead;
        stateHead_1 = stateHead;
        stateHead_2 = stateHead;

        for j = 1:size(stateHead_1, 1)
            if stateHead_2(j, 3) == 2
                stateHead_2(j, 2) = NaN;
            end
        end

        t = time_obs;
        plot(t, stateHead(:,2),'.-b','LineWidth',0.025);
        plot(t,stateHead_1(:,2),"Marker",".","Color","r");
        plot(t,stateHead_2(:,2),"Marker",".","Color",'b');
        tspan = min(t):1:max(t);
        filt = day(tspan)==1 & month(tspan)==1;
        tspan = tspan(filt);
        yspan = year(tspan);
        xticks(gca,tspan);
        xticklabels(gca,yspan);
        title(['ID' filename ' Two-State Solved Plot Compared to Single-State Solved Plot']);
        legend('Single State Model','Observed Groundwater Head','State 1','State 2');
        hold off

        % ADD Additional h_stars from model_6params.model.variables.h_stars
    end
end

if fourPlotsOneAxes == false
    for k = 2 %1:length(list)
        filename = list{k};

        % Load Two State model
        load(filename + "_model_6params_2.mat");
        time_steps = model_6params.simulationResults{1,1}.head(:,1);
        head_obs = model_6params.model.inputData.head(:,2);
        head_model1 = model_6params.model.variables.h_stars(:,2); %hstars(:,2)
        head_model2 = model_6params.model.variables.h_stars(:,3); %hstars(:,3)

        figure(k);
        plot(time_steps,head_obs,'k');
        hold on;

        dryIndex = model_6params.model.variables.StateHead;
        plot(time_steps, dryIndex(:,2),'.-b','LineWidth',0.025);

        for j = 1:size(dryIndex, 1)
            if dryIndex(j, 3) == 1
                dryIndex(j, 2) = NaN;
            end
        end
        plot(time_steps(~isnan(dryIndex(:,2))), head_model1(~isnan(dryIndex(:,2))),'.-','Color',[0.3010 0.7450 0.9330],'LineWidth',0.025);

        plot(time_steps,dryIndex(:,2),"Marker",".","Color","r");
        tspan = min(time_steps):1:max(time_steps);
        filt = day(tspan)==1 & month(tspan)==1;
        tspan = tspan(filt);
        yspan = year(tspan);
        xticks(gca,tspan);
        xticklabels(gca,yspan);
        title(['ID' filename ' Two-State Solved Plot Comparing States After Switch']);
        legend('Observed Groundwater Head','Wet State','Modelled Wet State','Dry State');
        hold off

        % ADD Additional h_stars from model_6params.model.variables.h_stars
    end
end