close all; clear; clc;

time_points = (726417:739386)';
delta_time = diff(time_points);

obshead = 250 + 10 .* randn(length(time_points),1);
obshead = [time_points,obshead];

h_star1 = 10 .* randn(length(time_points),1) + 249;
h_star1 = [time_points,h_star1];
h_star2 = 10 .* randn(length(time_points),1) + 251;
h_star2 = [time_points,h_star2];

alpha = log10(-log(eps(0)));

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

% Method 1 - PDF
sigma_n_1 = sqrt(mean( innov1.^2 ./ (1 - exp( -2 .* 10.^alpha .* delta_time ))));
sigma_n_2 = sqrt(mean( innov2.^2 ./ (1 - exp( -2 .* 10.^alpha .* delta_time ))));
objFn_1_PDF = pdf('Normal', resid1, 0, sigma_n_1);

objFn_2_PDF = pdf('Normal', resid2, 0, sigma_n_2);

% Method 5 - original objFn
original = false;
if original == true
    objFn_1_logPDF = normlike([0,sigma_n_1],resid1);
    objFn_1_sumlogPDF = -sum(log(pdf('Normal', resid1, 0, sigma_n_1)));
    N = size(resid1,1);
    objFn =  sum(exp(mean(log( 1- exp( -2.*10.^alpha .* delta_time) ))) ...
        ./(1- exp( -2.*10.^alpha .* delta_time )) .* innov1.^2);
    objFn = 0.5 * N * ( log(2*pi) + log(objFn/N)+1);
else
    sigma_n_1 = sqrt(mean( innov1(1).^2 ./ (1 - exp( -2 .* 10.^alpha .* delta_time(1) ))));
    objFn_1_logPDF = normlike([0,sigma_n_1],resid1(1));
    objFn_1_sumlogPDF = -sum(log(pdf('Normal', resid1(1), 0, sigma_n_1)));
    N = 1;
    objFn =  exp(mean(log( 1- exp( -2.*10.^alpha .* delta_time(1)) ))) ...
        ./(1- exp( -2.*10.^alpha .* delta_time(1) )) .* innov1(1).^2;
    objFn = 0.5 * N * ( log(2*pi) + log(objFn/N)+1);
end

test2 = (2.*pi() .* exp(1) * (exp(mean(log( 1- exp( -2.*10.^alpha .* delta_time) ))) ...
        ./(1- exp( -2.*10.^alpha .* delta_time )) .* innov1.^2)).^-0.5;

test_x = exp(mean(log(1-exp(-2.*10.^alpha .* delta_time))));
test = (2.* pi .* exp(1) .* innov1.^2 ./ (1 - exp(-2.*10.^alpha .* delta_time)) .* test_x).^-0.5;

% My findings
    % FINDING #1
    % When resid = innov, and we use the original objFn including the
    % sum(), the value it produces matches loglikelihood function that
    % matlab uses, specifically normlike with a mean of 0, and a stdev of
    % sigma_n_1 = sqrt(mean( innov1.^2 ./ (1 - exp( -2 .* 10.^alpha .* delta_time ))))
    % except it solves the negative loglikelihood, so the magnitude is the
    % same, the sign is different
    % So we're calculating the loglikelihood correctly, just not the log
    % probabilities correctly

    % FINDING #2
    % When resid = innov, and we use the original objFn including the
    % sum(), the value it produces matches loglikelihood function that
    % matlab uses, even if the only point of comparison is one datapoint.
    % However, in order for it to work, the datapoint(s) must also be used
    % to find sigma_n for the normlike

    % FINDING #3
    % We would expect that if innov =/= resid, the first value produced by
    % an objective function that considered historical correlation, would
    % match the first value produced by a PDF that does not consider
    % historical correlation, because innov(1) == resid(1) as it has no
    % history to consider