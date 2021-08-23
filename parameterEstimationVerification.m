clear 
close all
clc


 set(0,'DefaultAxesFontSize',18,'DefaultTextFontSize',24,...
     'DefaultAxesFontName','Helvetica',...
    'DefaultTextFontName','Helvetica',...
     'DefaultAxesFontWeight','bold','DefaultTextFontWeight','bold',...
     'DefaultLineLineWidth',2,'DefaultLineMarkerSize',12,...
     'DefaultFigureColor','w','DefaultFigureResize','on')

%  This Matlab script plots the cumulative and daily reported deaths 
%  obtained from JHU against mathematical model.

 %% US Data - from JHU 
% These are the number of deaths reported in US daily from March 1st to November 26th. 
deathsUSOrig = [1	0	5	1	4	1	2	3	4	1	6	5	10	8	7	12	27	35	59	73	101	99	110	188	241	326	412	518	633	592	696	1092	1180	1503	1393	1541	1568	1776	2570	2154	2219	2227	2124	1815	1938	2454	2603	2169	2093	1969	1931	2180	2542	2438	2466	2160	1703	1334	1465	2241	2527	2416	1877	1687	1103	1318	2299	2344	1924	1752	1542	900	1002	1614	1749	1774	1664	1272	755	1167	1526	1518	1211	1245	1134	616	557	664	1508	1183	1155	954	689	777	1049	994	1005	906	697	441	497	930	899	854	820	773	313	383	821	731	708	662	606	293	384	842	756	545	621	502	278	378	648	692	720	659	269	289	361	1207	863	998	821	693	461	368	920	972	950	933	862	453	524	1101	1220	1092	1130	888	483	1119	1369	1431	1216	1242	1109	411	539	1374	1376	1248	1238	1082	516	527	1068	1506	1064	1343	1030	575	446	1331	1331	1078	1111	984	455	444	1230	1228	1104	972	964	312	595	1065	1052	1067	970	771	410	276	446	1182	907	1213	710	386	424	1288	987	866	930	715	241	425	1023	1102	917	939	746	282	321	912	942	851	928	679	349	462	702	930	964	965	627	422	337	797	1001	818	905	745	428	467	934	1155	850	951	936	378	492	968	1007	967	1028	885	408	532	1571	1090	1150	1195	1008	473	688	1381	1413	1209	1158	1260	627	734	1692	1897	1987	1917	1429	920	913	2146	2297	1232];

% Population Data 
N = 328.2e6; % US population

%% For plotting purposes. 

% % Selectingdata range
begD = 1;  %March 1st
endD = numel(deathsUSOrig);

% Date formating. 
startDate = datetime(2020,3,1);
otherDates = startDate+caldays(1:endD-begD);
allDates = [startDate,otherDates];


%% Processing of data obtained from JHU. 

% obtaining 7-day moving average. 
deathsDailyData = movmean(deathsUSOrig,7); % 7-day moving average based on daily reported deaths. 
% Cumulative deaths based on moving average of daily reported deaths
deathsCumData = cumsum(deathsDailyData);

%% Simulation results

tspan = [0 400];
times = 1:1:400;

% Model parameters
%      N       Sm0     Eu0 Em0 Iu0 Im0 Au0 Am0 Ih0 Iicu0   R0      D0   beta1    beta2   beta3     beta4        beta5   Lambda1     lambda2      lambda3  theta1  theta2  tau1    tau2    tau3    tau4     
fixp = [N 3276578 0   0   169 0   440 877 0   0       2014    1    0.225    0.659   0.097     0.567112    0.10658 0.179      13.53       0.0388   0.1098  0.847   27.5    97.98   127.6   222.13];

%           Cbeta Ctheta    epsI    epsO    betaF   muF     F      etaA    etaH    Nu       Sigma   phi     gammaI  gammaA      gammaH      gammaICU    deltaI      deltaA      deltaH      deltaICU 
varpNom = [ 1     1         0.5     0.5     4.8     0.54    0.25   0.5     0.5     0.083    1/5.1   1/5     1/10    0.13978     1/8         1/10        0.015       0.0075      0.015       0.0225  ];
paraName = {'C_\beta','C_\theta','\epsilon_I','\epsilon_O','\beta_F','\mu_F','f','\etaA','\eta_H','\nu','\sigma','\phi','\gamma_I','\gamma_A','\gamma_H','\gamma_{ICU}','\delta_I','\delta_A','\delta_H','\delta_{ICU}'};

% Model Calculation 

% SEIR12EqnCalc.m provides all the model equations and how time dependent
% transmission coefficients are implemented. 

% IC12.m assins the initial conditions required for the simulation.
opt   = odeset('RelTol', 1e-7, 'AbsTol', 1.0e-9); % ode45 options
[x] = ode15s(@SEIR12EqnCalc,tspan,IC12(fixp),opt,varpNom,fixp);
xpp = deval(x,times);
deathsCumModel = xpp(12,:);
deathsDailyModel(1) = deathsCumModel(1);
for j = 2:numel(deathsCumModel)
    deathsDailyModel(j) = deathsCumModel(j)-deathsCumModel(j-1); 
end
allDatesM = [startDate,startDate+caldays(1:399)];
allDatesPlot = allDatesM(times);

%% Plotting 

% this provides a comparison of original JHU results and the 7-day moving
% average results. 
figure(1)
plot(allDates, deathsUSOrig(begD:endD), 'rs-')
hold on;
plot(allDates, deathsDailyData, 'ko')
legend('Daily Deaths - US','Daily Deaths (7-day Moving Average) - US','Location','southoutside')
xlabel('Time (Days)')
legend('boxoff')   
ylabel('Number of Deaths')
set(gcf, 'Position', [10 10 1200 675]);


% comparison of JHU moving average data and the model results. 
figure(2)
subplot(2,1,1);
plot(allDates, deathsDailyData, 'ko')
hold on;
plot(allDatesM(1:306), deathsDailyModel(1:306),'r-')
legend('Daily Deaths','Model','Location','northwest')
xlabel('Time (Days)')
ylabel('Daily Deaths')

subplot(2,1,2); 
plot(allDates, deathsCumData, 'ko')
hold on;
plot(allDatesM(1:306), deathsCumModel(1:306),'r-')
legend('Cum. Deaths','Model','Location','northwest')
xlabel('Time (Days)')
ylabel('Cumulative Deaths')
set(gcf, 'Position', [10 10 1000 750]);