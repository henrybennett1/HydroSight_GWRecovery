close all; clear; clc;

alpha = inf;

time_points = (726417:739386)';
delta_time = diff(time_points);

obshead = 250 + 10 .* randn(length(time_points),1);
obshead = [time_points,obshead];

h_star1 = 10 .* randn(length(time_points),1) + 249;
h_star1 = [time_points,h_star1];
h_star2 = 10 .* randn(length(time_points),1) + 251;
h_star2 = [time_points,h_star2];

t_filt = find( obshead(:,1) >=time_points(1)  ...
    & obshead(:,1) <= time_points(end) );
resid1= obshead(t_filt,2)  - h_star1(:,2);
resid2= obshead(t_filt,2)  - h_star2(:,2);
delta_time = [inf;delta_time];
% Calculate innovations using residuals from the deterministic components.
innov1 = resid1(2:end) - resid1(1:end-1).*exp( -10.^alpha .* delta_time(2:end) );
innov2 = resid2(2:end) - resid2(1:end-1).*exp( -10.^alpha .* delta_time(2:end) );
% at first timestep prior contributions is zero
innov1 = [resid1(1);innov1];
innov2 = [resid2(1);innov2];

% Requirements for objFn to be considered successful

% When alpha = inf || delta_time = inf, innov == resid
% When innov == resid, 
%       objFn == pdf(resid, 0, stdev)
%       Where stdev = sqrt(mean( innov1.^2 ./ (1 - exp( -2 .* 10.^alpha .* delta_time ))));

% For all values of alpha
%       sum(log(objFn)) == -original_objFn


% TESTING
% WHAT IS sigma_n, sigma_v and sigma_a, how are they different ?
