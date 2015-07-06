function TSDRcomparisonPlots(u,costs,uD,costsD,w,iFigure, params)
%TSDRcomparisonPlots compares policies with and without the demand charge
%in the TSDR problem.
%
% This function plots the net electrical load and costs of daily operation
% for the policies with and without the demand charge explicitly included.
%
% ---- Inputs ----
%
%   u               =[u(0),...,u(N-1)], the nu x N matrix of controls,
%                   calculated with regard to the demand charge.
%
%   costs           a struct containing the energy, DR, disutility,
%                   demand and terminal costs, calculated with regard to
%                   the demand charge.
%
%   uD              =[uNaive(0),...,uNaive(N-1)], the nu x N matrix of 
%                   controls, calculated without regard to the demand
%                   charge.
%
%   costsD          a struct containing the energy, DR, disutility,
%                   demand and terminal costs, calculated without regard to
%                   the demand charge.
%
%   w               =[w(0),...,w(N-1)]is the true disturbance.
%
%   iFigure         the integer index of the first figure to plot into.
%
% params is a struct containing:
%
%   N               the planning horizon.
%
%   dt              the problem time step.
%
%   kMain           the coefficient of performance of the main chiller.
%
%   pBase           the demand response baseline.


%
% Import the horizon and define a time vector to plot over.
%
  N = params.N;
  dt = params.dt;
  ks = 1 : N;
%
% Define the title font size.
%
  titleSize = 16;
%

%% power plots
% This cell plots the total power for both policies.

%
% Plot the power calculated without the demand charge.
%
  figure(iFigure); clf
  powerAxes1 = subplot(2,2,1);
  hold on
  stairs(ks, uD(1,:) + uD(3,:) + w(2,:), 'b')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2), ...
    'YLim',[0,50*ceil(max(uD(1,:) + uD(3,:) + w(2,:))/50)], ...
    'Box', 'on')
  title('Demand Charge Ignored')
  ylabel('Total Power (kW)')
  xlabel('Time (hours)')
%
% Plot the power calculated with the demand charge.
%
  powerAxes2 = subplot(2,2,2);
  hold on
  stairs(ks, u(1,:) + u(3,:) + w(2,:), 'r')
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2), ...
    'YLim',[0,50*ceil(max(uD(1,:) + uD(3,:) + w(2,:))/50)], ...
    'Box', 'on')
  title('Demand Charge Included')
  ylabel('Total Power (kW)')
  xlabel('Time (hours)')
%


%% cost plots 
% This cell plots the costs of both policies.

%
% Compute the net energy, curtailment, demand, disutility, and 
% noncompliance costs, calculated without the demand charge. Put them in 
% arrays with the terminal cost.
%
  netEnergyCostD = sum(costsD.energyCosts);
  netCurtailmentCostD = sum(costsD.DRcosts);
  demandCostD = costsD.demandCost;
  netDisutilityCostD = sum(costs.disutilityCosts);
  tankCostD = costsD.tankCost;
  costArrayD = [netEnergyCostD, netCurtailmentCostD, demandCostD, ...
               netDisutilityCostD, tankCostD];
%
% Compute the net energy, curtailment, demand, disutility, and 
% noncompliance costs, calculated with the demand charge. Put them in 
% arrays with the terminal cost.
%
  netEnergyCost = sum(costs.energyCosts);
  netCurtailmentCost = sum(costs.DRcosts);
  demandCost = costs.demandCost;
  netDisutilityCost = sum(costs.disutilityCosts);
  tankCost = costs.tankCost;
  costArray = [netEnergyCost, netCurtailmentCost, demandCost, ...
               netDisutilityCost, tankCost];
%
% Plot the cost array in a bar graph in the lower two subplots.
%
  barAxes = subplot(2,2,3:4);
  bar([costArrayD', costArray'])
%
% Get the upper subplot dimensions and align the lower subplot.
%
%   power1pos = get(powerAxes1,'Position');
%   power2pos = get(powerAxes2,'Position');
%   barPos = get(barAxes,'Position');
%   barPos(1) = power1pos(1);
%   barPos(3) = power1pos(3) + power2pos(3) + ...
%     power2pos(1) - (power1pos(1) + power1pos(3));
%   set(barAxes,'Position',barPos);
%
% Set the x- and y-axis bounds.
%
  yMin = roundn(1.1*netCurtailmentCostD,2);
  yMax = roundn(1.1*demandCostD,2);
  xMin = 0;
  xMax = length(costArrayD)+1;
  set(gca,'YLim', [yMin,yMax], 'XLim', [xMin,xMax]);
  title('Costs with and without Demand Charge in Objective Function')
%
% Set the x-tick labels.
%
  set(gca,'XTickLabel',{'Energy', 'Demand Response', 'Demand', ...
    'Under-cooling', 'Tank Depletion'});
%
% Add a y-axis label and a legend.
%
  ylabel('\$')
  legend('Without $g_d$ \quad','With $g_d$', 'Location', 'NorthEast')
%


%% end of function TSDRcomparisonPlots


end

