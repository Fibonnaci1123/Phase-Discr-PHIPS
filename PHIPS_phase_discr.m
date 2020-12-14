function [dropletflag_algorithm] = PHIPS_phase_discr(ASF);
% This matlab code determines the phase of a single particle based on its ASF
% measured by PHIPS
% input: ASF = the scattering intensities of the 20 nephelometer channels
% output: dropletflag, 1 = droplet, 0 = ice
% last updated: 14.12.2020
% copyright Fritz Waitz, KIT, IMK-AAF, fritzwaitz@kit.edu


f1 = calculate_f1(ASF);
f2 = calculate_f2(ASF);
f3 = calculate_f3(ASF);
f4 = calculate_f4(ASF);

f = [f1,f2,f3,f4];

[dropletflag_algorithm, confidence ,prob] = calculate_probability(f);
end

%% Divide by mie and integrate the difference
function f1 = calculate_f1(ASF);
angle = 18:8:170;

%++++++ 
% calculated Mie intensity for a droplet D = 100um, integrated over nephelometer FoV
I_int = [1.38832870096382e-05,8.66487003343236e-06,5.29017936106249e-06,2.86841609162495e-06,1.66423607932730e-06,8.48069104065888e-07,3.38763124347349e-07,1.24602153109646e-07,5.64338044300031e-08,3.56481856823524e-08,3.71664129976456e-08,3.27267735708208e-08,4.02623240663777e-08,6.60656796279482e-08,1.12817985291965e-07,5.64320582926844e-07,4.54212580302176e-07,2.61223649393121e-07,1.87147367956421e-07,1.57711848996302e-07];
%++++++    

% divide by mie
mie_ratio = ASF./I_int;

% delete 90 and NaN's (e.g. 34 in ACLOUD)
idx = find(isnan(mie_ratio));
angle(idx) =  [];
mie_ratio(:,idx) = [];

if length(angle) <= 5
    f1 = NaN;
else
    % normalize by mean of SPF
    % log makes it so that /2 has the same weight as x2
    % also, due to the log, it is now centered around 0 instead of 1
    mie_ratio = log(mie_ratio ./ median(mie_ratio,2));
    % integrate
    f1 = log ( trapz(angle, abs(mie_ratio),2) );
end
end

%% Downward Slope
function f2 = calculate_f2(ASF);
%%
angle = 18:8:170;
angle_down_start = 42;
angle_down_end = 74;

fit_angle = (angle_down_start -10)/8:(angle_down_end -10)/8;
fit_angle(isnan(ASF(fit_angle))) =  [];% remove those NaNs from index

if size(fit_angle,2) < 2 %skip this particle if only 1 angle is =/= NaN
    f2 = NaN;
else
    p = polyfit(angle(1,fit_angle),log(ASF(fit_angle)),1);
    pf = zeros(1,size(fit_angle,2));
    
    for j=1:size(fit_angle,2)
        pf(j) = exp ( angle(1,fit_angle(j)).* p(1) +  p(2));
    end
    f2  = (log(pf(end)) - log(pf(1))) ./ (angle(1,fit_angle(end))-angle(1,fit_angle(1))) * 100; %*100 so that the numbers are nicer
end
end

%% Upward Slope
function f3 = calculate_f3(ASF);
angle = 18:8:170;

angle_up_min_start = 114;
angle_up_min_end = 130;
angle_up_max_start = 130;
angle_up_max_end = 146;
angle_up_min = (angle_up_min_start -10)/8:(angle_up_min_end -10)/8; %define angle interval
angle_up_max = (angle_up_max_start -10)/8:(angle_up_max_end -10)/8; %define angle interval

angle_up_min(isnan(ASF(angle_up_min))) =  [];% remove those NaNs from index
angle_up_max(isnan(ASF(angle_up_max))) =  [];% remove those NaNs from index

if size(angle_up_min,2) < 1 || size(angle_up_max,2) < 1 %skip this particle if there is no max or min angle       
    f3 = NaN;
else       
    [u_min, position_min] = min (ASF(angle_up_min)); %find minimum
    [u_max, position_max] = max (ASF(angle_up_max));
    fit_angle = angle_up_min(1) + position_min - 1 : angle_up_max(1) + position_max - 1; %interval from min to max
   
    if fit_angle(1) == fit_angle(end)+1 || fit_angle(1) == fit_angle(end)-1 %if max = min +/- 1
        fit_angle = angle_up_min(1):angle_up_max(end);       % then use whole interval
    end  

    fit_angle(isnan(ASF(fit_angle))) =  [];% remove those NaNs from index

    p = polyfit(angle(1,fit_angle),log(ASF(fit_angle)),1); 
    pf = zeros(1,size(fit_angle,2)); 
    for j=1:size(fit_angle,2)
        pf(j) = exp (angle(1,fit_angle(j)).* p(1) +  p(2));
    end

    f3 = (log(pf(end)) - log(pf(1))) ./ (angle(1,fit_angle(end))-angle(1,fit_angle(1))) * 100; %*100 so that the numbers are nicer

    if fit_angle(1) > fit_angle(end) %if max_angle < min_angle, fix the sign 
        f3 = - f3;
    end

    if (fit_angle(1) == fit_angle(end)+1 || fit_angle(1) == fit_angle(end)-1) && f3 < 2 %up_ice_cond = 2 %if there is no min/max or if they are +/-1 close to each other
        f3 = NaN; 
    end
end

end

%% ratio around the minimum at 90
function f4 = calculate_f4(ASF);
ang = [8,9,11,12];
a = mean(ASF(ang));
bang = [6,7,13,14];
b = mean(ASF(bang));
f4 = log ( a/b );
end

%% calculate probabilities and dropletflag
function [dropletflag_algorithm ,confidence,prob] = calculate_probability(f);
% f = [f1,f2,f3,f4);
weight = [1,1,1,1];
% use the calibration from ACLOUD for all campaigns
prob_fit_coeff =[   98.5671   10.8862    4.6407    0.1804
                   99.3590    2.5199   -7.3120    2.0520
                  98.0423   -1.0686    5.7315    3.1399
                 99.2689   10.7820   -0.9897    2.1943
                ];

% % Calculate ice/droplet probability
for i = 1:4
    a = prob_fit_coeff(i,1);
    b = prob_fit_coeff(i,2);
    c = prob_fit_coeff(i,3);
    d = prob_fit_coeff(i,4);

    prob(i) = (a-d)./ (1 + exp(-b.*(f(i)-c))) + d ;

end

weighted_prob = weight .* prob;
% check which f (and hence p) are NaN, 1 if ok, 0 if NaN
% delete the weighting from NaN entries
weight_array = abs(isnan(f)-1) .* weight; % if a slope is NaN, then for this particle, the weight for this slope is zero

confidence =  100 - nansum(weighted_prob,2)./nansum(weight_array,2); %confidence that a particle is ice, 100 = ice, 0 = droplet, 50 = not sure

dropletflag_algorithm(confidence>=50) = 1;
dropletflag_algorithm(confidence<50) = 0;
dropletflag_algorithm(isnan(confidence)) = NaN;
end
