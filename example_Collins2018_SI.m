% Collins2018 Supplementary Information Fig S8
% this is to check how they got their background-subtracted heatrates
% and to compare it with heatrates including the background buffer.
%
% Note that values are eye-balled from Fig S8 using a drawn horizontal line
% as aid, and therefore involve some uncertainty (\pm 0.1, maybe 0.2).
% Sample C NP+buffer values are especially uncertain due to few y-axis
% ticks.

t = 300; % [s]

%% Sample A: pH 5
% NP+buffer
dT_npb = [32.9-20.8, 32-20, 30.9-19.3];
dT_npb = mean(dT_npb);
% only buffer 
dT_b = [26-20.6, 25-19.3, 24-19];
dT_b = mean(dT_b);
% background+NPs heatrate
dT_A_all = dT_npb/t
% background-subtracted heatrate
dT_A = (dT_npb-dT_b)/t

%% Sample B: pH 7
% NP+buffer
dT_npb = [34-22, 32.9-22, 31.4-22];
dT_npb = mean(dT_npb);
% only buffer 
dT_b = [24.6-22.25, 24.5-21.9, 24.3-21.9];
dT_b = mean(dT_b);
% background+NPs heatrate
dT_B_all = dT_npb/t
% background-subtracted heatrate
dT_B = (dT_npb-dT_b)/t

%% Sample C: pH 9
% NP+buffer
dT_npb = [32-19.5, 31.5-19.4, 30-18];
dT_npb = mean(dT_npb);
% only buffer 
dT_b = [21.5-18.2, 20-17.3, 20-17.4];
dT_b = mean(dT_b);
% background+NPs heatrate
dT_C_all = dT_npb/t
% background-subtracted heatrate
dT_C = (dT_npb-dT_b)/t