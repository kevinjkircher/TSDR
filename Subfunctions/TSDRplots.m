function TSDRplots(x,u,w,costs,iFigure,params,legLabel)
%TSDRplots draws time series and histograms for the TSDR problem.
%
% This function plots time series of the states and controls in one figure.
% It plots time series of the total cooling power alongside the prices of
% energy, demand and curtailment in another figure. In a third figure, it
% plots the stage-wise and total costs.
%
% ---- Inputs ----
%
%   x               =[x(0), ..., x(N)] is the nx x N+1 matrix of state
%                   vectors.
%
%   u               =[u(0), ..., u(N-1)] is the nu x N matrix of control
%                   vectors.
%
%   w               =[w(0), ..., w(N-1)] is the nw x N matrix of
%                   disturbance vectors.
%
%   iFigure         the integer index of the first figure to plot into.
%
%   legLabel        the legend entry (typically oracle, baseline or MPC), a
%                   string.
%
% costs is a struct containing:
%
%   energyCosts         A 1 x nTimes vector of the energy costs at stages
%                       0, 1, ..., nTimes-1.
%
%   DRcosts             A 1 x nTimes vector of the DR costs at stages 
%                       0, 1, ..., nTimes-1.
%
%   disutilityCosts     A 1 x nTimes vector of the disutility costs at 
%                       stages 0, 1, ..., nTimes-1.
%
%   tankCost            The scalar cost of tank depletion.
%
%   demandCost          The scalar peak demand cost.
%
%   totalCost           The scalar total cost.
%

%
% Import the horizon.
%
  N = params.N;
%

%% state and control time series 
% This cell plots the states and controls.

%
% Plot x1, the tank charge state.
%
  ks = 0 : 1 : N;
  figure(iFigure); clf
  subplot(5,1,1) 
  stairs(ks,x(1,:))
  title('Tank Charge State')
  ylabel('$x_1$ (kWh$_{th})$')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%
% Plot x2, the backlogged cooling load.
%
  subplot(5,1,2) 
  stairs(ks,x(2,:))
  title('Cooling Load Deficit')
  ylabel('$x_2$ (kWh$_{th}$)')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%
% Plot u1, the power consumed by the ice chiller.
%
  subplot(5,1,3) 
  stairs(ks(1:end-1),u(1,:))
  title('Power to Ice Chiller')
  ylabel('$u_1$ (kW)')  
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%
% Plot u2, the cooling provided by melting ice.
%
  subplot(5,1,4) 
  stairs(ks(1:end-1),u(2,:))
  title('Cooling from Ice Melt')
  ylabel('$u_2$ (kW$_{th}$)')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%
% Plot u3, the power consumed by the main chiller.
%
  subplot(5,1,5) 
  stairs(ks(1:end-1),u(3,:))
  title('Power to Main Chiller')
  ylabel('$u_3$ (kW)')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
  xlabel('Time (hours)')
%


%% electricity and cost time series; cost bar graph
% This cell plots the electricity consumption and its baseline, as well as
% the stage costs and a bar graph of the total cost.

%
% Plot the total electric power, its baseline, and the power that would be
% consumed if all cooling load were met directly with the main chiller
% (called pNaive).
%
  pNaive = w(1,:)./params.kMain + w(2,:);
  figure(iFigure+1); clf
  subplot(2,1,1)
  hold on
  stairs(ks(1:end-1), u(1,:) + u(3,:) + w(2,:), 'k.-')
  stairs(ks(1:end-1), params.pBase, 'r')
  stairs(ks(1:end-1), pNaive, 'c')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2), ...
    'YLim',[0,50*ceil(max(pNaive)/50)])
  title('Total Power Consumption')
  ylabel('kW')
  legend(legLabel, 'Baseline', 'No Storage \quad', ...
    'Location', 'SouthEast', 'Orientation', 'Horizontal')
  box on
% %
% % Plot the energy and DR prices.
% %
%   subplot(3,1,3); hold on
%   [edAx,h1,h2] = plotyy(ks(1:end-1), params.ce, ks(1:end-1), params.cDR, ...
%    'stairs', 'stairs');
%   title('Energy and Curtailment Prices')
%   ylabel(edAx(1),'$c_e(k)$ (\$/kWh)')
%   set(edAx(1),'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%   ylabel(edAx(2),'$c_{dr}\ (k)$ (\$/kWh)')
%   set(edAx(2),'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
% %
% % Plot the disutility price.
% %
%   subplot(4,1,4)
%   stairs(params.cu, 'r')
%   title('Price of Under- or Over-cooling')
%   ylabel('\$/kWh$_{th}^2$')
%   set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
%   xlabel('Time (hours)')
%
%
% Plot the stage costs.
%
%  figure(iFigure+2); clf
  subplot(2,1,2)
  barsToPlot = [costs.energyCosts', costs.DRcosts', costs.disutilityCosts'];
  bar(0:N-1, barsToPlot, 'stack')
  title('Costs')
  xlabel('Time (hours)')
  ylabel('\$')
  legend('Energy', 'Demand Response\quad', 'Under-cooling', ...
    'Location', 'SouthWest')
  set(gca, 'XLim', [0,N], 'XTick',0:12:N,'XTickLabel',0:6:round(N/2))
% %
% % Compute the net energy, curtailment, demand, disutility, and 
% % noncompliance costs. Put them in an array with the tank cost.
% %
%   netEnergyCost = sum(costs.energyCosts);
%   netCurtailmentCost = sum(costs.DRcosts);
%   demandCost = costs.demandCost;
%   netDisutilityCost = sum(costs.disutilityCosts);
%   costArray = [netEnergyCost, netCurtailmentCost, demandCost, ...
%                netDisutilityCost, costs.tankCost];
% %
% % Plot the cost array with text labels on the x axis.
% %
%   subplot(3,1,3)
%   bar(costArray)
%   set(gca,'XTickLabel',{'Energy', 'Demand Response', 'Demand', ...
%                 'Under-cooling', 'Tank Depletion'});
%   title('Cost of Daily Operation')
%   ylabel('\$')
% %


%% end of function TSDRplots


end

