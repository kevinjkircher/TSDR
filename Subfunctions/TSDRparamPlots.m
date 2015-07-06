function TSDRparamPlots(iFigure, params)
%TSDRparamPlots plots the time-varying parameters in the TSDR problem.
%
% This function plots time series of the temperature, chiller coefficients
% of performance, energy costs, demand costs, costs of unmet (or over-met)
% cooling load, curtailment costs, and the expected loads and their standard
% deviations.
%
% ---- Inputs ----
%
%   iFigure         the integer index of the first figure to plot into.
%
% params is a struct containing:
%
%   temperatures    the temperature at each stage, in Farenheit.
%
%   kMain           the coefficients of performance of the main chiller.
%
%   kIce            the coefficients of performance of the ice chiller.
%
%   ce              the energy costs, in $/kWh.
%
%   cd              the demand costs, in $/kW.
%
%   cu              the costs of unmet cooling load, in $/kWhth^2.
%
%   cDR             the demand response curtailment costs, in $/kWh.
%
%   wbar            =[w1bar; w2bar], the expected cooling loads and firm
%                   loads, in kWth and kW, respectively.
%
%   Q               =[Q(0),...,Q(N-1)], the load covariance matrices.
%


%
% Import the horizon and define a time vector to plot over.
%
  N = params.N;
  ks = 1 : N;
%
% Define the title font size.
%
  titleSize = 16;
%

%% temperature, COP and load time series 
% This cell plots the temperatures and coefficients of performance in one
% subplot, and the expected loads in another.

%
% Plot the temperature and coefficients of performance in the top subplot.
% Adjust the line and axis colors.
%
  figure(iFigure); clf
  subplot(2,1,1); hold on
  [tempAx,h1,h2] = plotyy(ks, (params.temperatures-32)*5/9, ks, params.kMain, ...
    'stairs', 'stairs');
  title('Temperature and Coefficients of Performance', 'FontSize', titleSize)
  set(h1,'Color','Red')
  set(tempAx(1),'YColor','Black')
  set(h2,'Color','Blue')
  set(tempAx(2),'YColor','Black')
%
% Add the ice chiller COP.
%
  axes(tempAx(2)); hold on
  stairs(ks, params.kIce, 'c')
%
% Add axis labels and adjust axis bounds.
%
  ylabel(tempAx(1),'Temperature ($^\circ$C)')
  set(tempAx(1),'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
  ylabel(tempAx(2),'COP')
  set(tempAx(2),'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2),...
    'YLim',[0,4],'YTick',0:1:4)
%
% Add a legend.
%
  tempLeg = legend('Main Chiller', 'Ice Chiller', 'Temperature\quad', ...
    'Location', 'SouthEast');
  set(tempLeg, 'Color', 'none');
%
% Plot the loads and baseline in the bottom subplot. Start by defining
% vectors of the standard deviations, pulled from Q(k).
%
  subplot(2,1,2); hold on
  [loadAx,h1,h2] = plotyy(ks, params.wbar(1,:), ...
    ks, params.wbar(2,:), 'stairs', 'stairs');
  title('Expected Loads and 95\% Confidence Intervals', 'FontSize', titleSize)
  set(h1,'Color','Blue')
  set(loadAx(1),'YColor','Blue','YLim',[0,500],'YTick',0:100:500,...
    'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
  set(h2,'Color','Red')
  set(loadAx(2),'Ycolor','Red','YLim',[0,500],'YTick',0:100:500,...
    'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%
% Add axis labels and adjust axis bounds.
%
  ylabel(loadAx(1),'Cooling Load (kW$_{th}$)')
  ylabel(loadAx(2),'Other Electric Loads (kW)')
  xlabel('Time (hours)')
%
% Add error bars.
%
  stdev1 = squeeze(sqrt(params.Q(1,1,:)))';
  stdev2 = squeeze(sqrt(params.Q(2,2,:)))';
  axes(loadAx(1)); hold on
  plot(ks, params.wbar(1,:) + 2*stdev1, 'b+', ...
    ks, params.wbar(1,:) - 2*stdev1, 'b+');
  axes(loadAx(2)); hold on
  plot(ks, params.wbar(2,:) + 2*stdev2, 'r+', ...
    ks, params.wbar(2,:) - 2*stdev2, 'r+');
%


%% price plots 
% This cell plots the energy prices, demand prices, prices of mis-met 
% cooling load, and DR prices.

%
% Get the off-peak, shoulder and peak demand prices and build them into the
% appropriate length arrays for plotting.
%
  offPeakVec = params.cdOffPeak*ones(1,N);
  shoulderVec = zeros(1,N);
  shoulderVec(params.shoulderStart:params.shoulderEnd) = ...
    params.cdShoulder;
  peakVec = zeros(1,N);
  peakVec(params.peakStart:params.peakEnd) = ...
    params.cdPeak;
  
%
% In the first subplot, plot the energy and DR prices.
%
  figure(iFigure+1); clf
  subplot(3,1,1); hold on
  [edAx,h1,h2] = plotyy(ks, params.ce, ks, params.cDR, ...
   'stairs', 'stairs');
  title('Energy and Demand Response Prices', 'FontSize', titleSize)
  ylabel(edAx(1),'$c_e(k)$ (\$/kWh)')
  set(edAx(1),'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
  ylabel(edAx(2),'$c_{dr}\ (k)$ (\$/kWh)')
  set(edAx(2),'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%
% In the second subplot, plot the demand prices.
%
  subplot(3,1,2); hold on
  stairs(ks, offPeakVec, 'g')
  stairs(ks, shoulderVec, 'm')
  stairs(ks, peakVec, 'r')
  title('Demand Prices', 'FontSize', titleSize)
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2),...
    'YLim',[0,20],'YTick',0:5:20,'Box','on')
  ylabel('$c_d(\mathcal{T}_i)$ (\$/kW)')
  legend('$c_d(\mathcal{T}_1)$','$c_d(\mathcal{T}_2)$',...
    '$c_d(\mathcal{T}_3$)','Location','SouthWest')
%
% In the third subplot, plot the prices of unmet or overmet cooling load.
%
  subplot(3,1,3)
  stairs(ks, params.cu)
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
  title('Prices of Under- or Over-cooling', 'FontSize', titleSize)
  xlabel('Time (hours)')
  ylabel('$c_u(k)$ (\$/kWh$^2$)')
%


%% disutility vs. time and deficit
% This cell draws two plots in one figure. The first plot is a 3D graph of 
% the disutility cost as a function of time and cooling deficit. The second
% plot is a slice of the first plot during the occupied hours, shown with
% several values of alpha.

%
% Define the deficit values and alpha values to plot over.
%
  nDeficits = 100;
  maxDeficit = 500;
  x2 = linspace(-maxDeficit,maxDeficit,nDeficits);
  alphas = [0,500,1000];
  nAlphas = length(alphas);
  times = 1:N;
%
% Define the values of the disutility function to plot.
%
  guValues = zeros(nDeficits,N,nAlphas);
  for iAlphas = 1 : nAlphas
    params.alpha = alphas(iAlphas);
    params.gu = @(k,xk) params.cu(k)*(xk(2)^2 + ...
      params.alpha*max(0,xk(2)));
    for iTimes = 1 : N
      for iDeficits = 1 : nDeficits
        guValues(iDeficits,iTimes,iAlphas) = params.gu(iTimes,...
          [0;x2(iDeficits);0]);
      end
    end
  end
%
% Draw the 3D plot of disutility vs. time and deficit.
%
  figure(iFigure+2); clf
  subplot(3,1,1:2)
  mesh(times,x2,guValues(:,:,nAlphas))
  view(45,30)
%  set(gca,'ZLim',[0,200]);
  set(gca,'XLim',[0,48],'XTick',0:8:48,'XTickLabel',0:4:24)
  xlabel('Time (hours)')%,'Interpreter','tex')
  ylabel('$x_2$')%,'Interpreter','tex')
  zlabel('$g_u$')%,'Interpreter','tex')
%
% Define the time value to plot and draw a slice of the 3D plot for each
% alpha value.
%
  subplot(3,1,3); hold on
  kPlot = 20;
  plot(x2,guValues(:,kPlot,1),'g',x2,guValues(:,kPlot,2),'b',...
    x2,guValues(:,kPlot,3),'r')
  xlabel('$x_2$')%,'Interpreter','tex')
  ylabel('$g_u$')%,'Interpreter','tex')
%  title('Cost of Unmet Cooling Load at 10 AM','Interpreter','tex')
  legStrs = cell(1,3);
  for iAlphas = 1 : nAlphas
    legStrs{iAlphas} = sprintf('$\\alpha$ = %i \\quad',alphas(iAlphas));
  end
  legend(legStrs, 'Location', 'NorthWest')
  set(gca,'Box','on');
%


%% end of function TSDRparamPlots


end

