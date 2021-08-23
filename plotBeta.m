clear 
close all
clc


 set(0,'DefaultAxesFontSize',18,'DefaultTextFontSize',24,...
     'DefaultAxesFontName','Helvetica',...
    'DefaultTextFontName','Helvetica',...
     'DefaultAxesFontWeight','bold','DefaultTextFontWeight','bold',...
     'DefaultLineLineWidth',2,'DefaultLineMarkerSize',12,...
     'DefaultFigureColor','w','DefaultFigureResize','on')
 
 %% Input parameters. 
 
%        beta1  beta2   beta3   beta4   beta5   tau     tau2    tau3    tau4    lambda1 lambda2
param = [0.225  0.659   0.097   0.5671  0.10658 27.5    97.98   127.6   222.1   0.179   13.53 ];

t = 1:1:272;
startDate = datetime(2020,3,1);
otherDates = startDate+caldays(1:numel(t)-1);
allDates = [startDate,otherDates];


for i =1:numel(t)
    betaVal(i) = betaCalc(param,t(i));
end

figure(1)
plot(allDates, betaVal, 'k-')
ylim([0,1])
xlabel('Time (Days)')  
ylabel('Transmission Coefficient')
set(gcf, 'Position',[10 10 1000 500]);
saveas(gcf,sprintf('beta'),'epsc');

function beta = betaCalc(param,t)
    beta1 = param(1);
    beta2 = param(2);
    beta3 = param(3);
    beta4 = param(4);
    beta5 = param(5);
    tau     = param(6);
    tau2    = param(7);
    tau3    = param(8);
    tau4    = param(9);
    lambda1 = param(10);
    lambda2 = param(11);
    
    
    if t <= tau
        beta  = beta1+beta2;
    elseif (t > tau)&&(t <= tau2)
        beta = beta1 + beta2*exp(-lambda1*(t-tau));
    elseif (t > tau2)&&(t <= tau3)
         beta = beta1 + beta2*exp(-lambda1*(tau2-tau))+beta3*(1-exp(-lambda2*(t-tau2)));
    elseif(t > tau3)&&(t <= tau4)
         beta = beta1 + beta2*exp(-lambda1*(tau2-tau))+beta3*(1-exp(-lambda2*(tau3-tau2)))+beta4*exp(-lambda1*(t-tau3)); 
    else
         beta = beta1 + beta2*exp(-lambda1*(tau2-tau))+beta3*(1-exp(-lambda2*(tau3-tau2)))+beta4*exp(-lambda1*(tau4-tau3))+beta5*(1-exp(-lambda2*(t-tau4))); 
    end
end

