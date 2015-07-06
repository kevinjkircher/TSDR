%% ---- MPC OF THERMAL STORAGE FOR DEMAND RESPONSE ---- 
%
%               Kevin Kircher, Cornell MAE
%               June 19, 2014
%
% If you find this script useful, please cite the companion paper:
%
%	K.J. Kircher and K.M. Zhang. "Model Predictive Control of
%	Thermal Storage for Demand Response." In American Control Conference, 
%	Proceedings of the 33rd (2015).
%
% This script simulates a building. The building can make ice, store it,
% and use it later to provide cooling. The building operator pays a dynamic
% energy price and occasionally receives a load-curtailment price signal.  
% The cooling is scheduled using model predictive control.
%
% Time-dependent vectors are stored as columns of a storage matrix. 
% Each storage matrix is organized such that the jth column holds the 
% vector at stage j. Columns index time.
%
% Time-dependent matrices are stored as 'pages' in a three-dimensional
% 'droor' (a 3D array). Each 3D storage array is organized such that the
% jth page holds the matrix at stage j. Pages index time.
%
% SI units are used throughout. Units of kW and kWh refer to electric
% power, while units of kWc and kWhc refer to cooling load. The two are
% related by the coefficient of performance of whatever system is
% providing the cooling:
%
%       COP = (heat removed from space) / (work done by cooling system)
%
% Time is in hours. Money is in dollars.
%
% The states of the system are
%
%       x1, the amount of thermal energy stored in the ice tank (kWhc)
%       x2, the backlogged cooling load (kWhc)
%
% The controls are
%
%       u1, the power allocated to ice-making (kW)
%       u2, the cooling load met by melting ice (kWhc)
%       u3, the power allocated to the main chiller (kW)
%
% The disturbance is
%
%       w1, the building's cooling load (kWhc)
%       w2, the building's firm load (kWh)
%
% This script assumes the default text interpreter for figures is LaTeX. If
% not, axis labels may give parser warnings. To interpret with LaTeX by
% default, type set(0,'DefaultTextInterpreter','LaTeX').
%
% The Subfunctions folder needs to be on to your Matlab path.

%% initialization of fixed parameters 
% This cell initializes the physical and economic parameters that are fixed
% over all simulation days.

%
% Generate a struct TSDRparams containing the fixed parameters.
%
  run fixedTSDRparams
%

%% declaration of simulation day and data import 
% This cell declares the year, month and day to be simulated.

%
% Mount the input data struct.
%
  load TSDRdata
  hourlyData = TESdata.hourlyData;
  realTimeData = TESdata.RTdata;
  yearsStudied = {'2009', '2010', '2011', '2012', '2013'};
%
% Choose the year, month and day to simulate and define the initial
% simulation hour.
%
  yearStr = '2013';
  iYear = find(ismember(yearsStudied, yearStr));
  iMonth = 7;
  iDay = 18;
%
% Find the starting and ending indices of the simulation day within the
% hourly data.
%
  startHour = sprintf('%s/%s/%s 0:00', num2str(iMonth), num2str(iDay), ...
    yearStr(end-1:end));
  endHour = sprintf('%s/%s/%s 23:00', num2str(iMonth), num2str(iDay), ...
    yearStr(end-1:end));
  hourlyStartIndex = find(strcmp(hourlyData{iYear,1}, startHour));
  hourlyEndIndex = find(strcmp(hourlyData{iYear,1}, endHour));
%
% Find the starting and ending indices of the simulation day within the
% real-time data.
%
  RTstartStr = sprintf('%s/%s/%s', num2str(iMonth), num2str(iDay), ...
    num2str(yearsStudied{iYear}(end-1:end)));
  RTstartIndex = find(strcmp(realTimeData{iYear,1}, RTstartStr), 1);
  RTendIndex = find(...
    strcmp(realTimeData{iYear,1}, RTstartStr), 1, 'last');
%
% Check: print the results.
%
  fprintf('Start (hourly data): %s\nEnd (hourly data): %s\n', ...
    hourlyData{iYear,1}{hourlyStartIndex}, ...
    hourlyData{iYear,1}{hourlyEndIndex});
  fprintf('Start (real-time data): %s %s\nEnd (real-time data): %s %s\n', ...
    realTimeData{iYear,1}{RTstartIndex}, ...
    realTimeData{iYear,5}{RTstartIndex}, ...
    realTimeData{iYear,1}{RTendIndex}, ...
    realTimeData{iYear,5}{RTendIndex})
%
% Clean up the workspace.
%
  clear TESdata
%

%% definition of tunable parameters 
% This cell defines the parameters in the problem that can be chosen by the
% controller.

%
% Define the tank depletion cost, ct ($/kWhth), such that a terminal 
% tank charge state x1(N) of 1 kWhth less than x1(0) incurs a cost  of 
% ct dollars.
%
  TSDRparams.ct = 3*mean(hourlyData{iYear,9} ...
    (hourlyStartIndex:hourlyStartIndex+8)/1000);
%
% Define the constant by which the ratio of occupancy to thermal mass will
% be scaled, in order to compute the prices of unmet load, cu(k).
%
  disutilityConstant = 0.01;%0.0014;
%
% Define alpha, the asymmetry tuning parameter in the stage cost of under-
% or over-cooling.
%
  TSDRparams.alpha = 1e3;
%
% Define the MPC horizon.
%
  H = N;
%
% Define pBar, the vector of peak demands (kW) in each interval T_i so far
% during the month.
%
  TSDRparams.pBar = zeros(3,N);
  TSDRparams.pBar(:,1) = [200;200;200];
%
% Define the initial tank charge state.
%
  TSDRparams.x0 = zeros(5,1);
  initialChargeFraction = 0.2;
  TSDRparams.x0(1) = initialChargeFraction*TSDRparams.x1max;
%
% Define the controls from stage -1, for use in ramping constraints at
% stage 0.
%
  TSDRparams.x0(3:5) = [0;0;0];
%

%% calculation of deterministic, day-dependent parameters 
% This cell uses the day's data to calculate the values of the parameters
% that are deterministic, but that depend on the day being simulated.

%
% Modify the struct TSDRparams to include the deterministic, day-dependent
% parameters.
%
  run dailyTSDRparams
%

%% baseline calculation 
% This cell solves the open-loop optimal control problem with the nominal
% cooling load and fixed load (a certainty-equivalence assumption) and with
% no DR participation. This generates the baseline cooling values used
% in the oracle and MPC policies.

%
% Call the OLOC solver with nominal disturbances and no DR incentives.
%
  wBase = TSDRparams.wbar;
  TSDRparams.cDR = zeros(1,N);
  TSDRparams.gDR = @(k,uk,wk) 0;
  [xBase, uBase, baseCosts] = TSDROLOC(wBase, TSDRparams);
%
% Store the baseline in the TSDRparams struct.
%
  TSDRparams.pBase = uBase(1,:) + uBase(3,:) + wBase(2,:);
%
% Plot the results.
%
  TSDRplots(xBase,uBase,wBase,baseCosts,1,TSDRparams,'Baseline')
%

%% demand response prices and stage cost redefinition 
% This cell defines the prices for load curtailment in ConEd's Commercial
% System Relief Program, under the Voluntary Participation Option. This
% program pays customers $3/kWh for reducing load in response to a DR call,
% which is issued at least 21 hours before an event.

%
% Define the DR payments. Assume it's a DR day, with the event falling
% between 2 and 6 PM (as it would in the Chelsea area of Manhattan).
%
  TSDRparams.cDR = zeros(1,N);
  if N == 24
    DRstart = 15;
    DRend = 19;
    TSDRparams.cDR(DRstart:DRend) = 3;
  end
  if N == 48
    DRstart = 29;
    DRend = 36;
    TSDRparams.cDR(DRstart:DRend) = 3;
  end
%
% Define the DR cost function.
%
  TSDRparams.delta = @(k,uk,wk) TSDRparams.pBase(k) - ...
    (uk(1) + uk(3) + wk(2));
  TSDRparams.gDR = @(k,uk,wk) -TSDRparams.cDR(k)*...
    TSDRparams.delta(k,uk,wk)*TSDRparams.dt;
%
% Redefine the net stage cost.
%
  TSDRparams.gk = @(k,xk,uk,wk) TSDRparams.ge(k,uk,wk) + ...
                 + TSDRparams.gu(k,xk) + TSDRparams.gDR(k,uk,wk);
%

%% disturbance generation 
% This cell generates the cooling loads and firm loads, which are modeled
% as jointly normal and independent from stage to stage.

%
% Build the disturbance vector at each stage using the appropriate mean and
% covariance.
%
  rng(1)
  w = zeros(2,N);
  for iTimes = 1 : N
    w(:,iTimes) = mvnrnd(TSDRparams.wbar(:,iTimes)',TSDRparams.Q(:,:,iTimes))';
  end
%

%% parameter plots 
% This cell plots the physical and economic parameters.

%
% Pass the figure index and parameter struct to the parameter plotting
% function.
%
  TSDRparamPlots(8,TSDRparams)
%

%% oracle policy 
% This cell solves the open-loop optimal control problem with the true
% disturbance values known perfectly at stage 0 (the oracle information
% pattern). No causal policy can achieve a lower cost than the oracle.

%
% Call the OLOC solver with the true disturbance values.
%
  tic; [xOracle, uOracle, oracleCosts] = TSDROLOC(w, TSDRparams);
  oracleRunTime = toc;
%
% Plot the results.
%
  TSDRplots(xOracle,uOracle,w,oracleCosts,3,TSDRparams,'Oracle')
%

%% MPC policy 
% This cell implements the model-predictive control algorithm.

%
% Define the nominal disturbance values.
%
  wnom = TSDRparams.wbar;
%
% Call the MPC function. This function call is the meat of the script.
% It solves a convex, constrained optimization problems at each stage 
% k in {0, ..., N-1}. The kth problem has (N-k)*nu decision variables. With
% H=N=48, the total number of variables decided is 3528 and the runtime is
% about four minutes on a 2.8 GHz processor.
%
  tic;
  [x, u, pBar, costs] = TSDRMPC(H, w, wnom, TSDRparams);
  MPCtime = toc;
%
% Plot the results.
%
  TSDRplots(x,u,w,costs,1,TSDRparams,'MPC')
%

%% MPC policy ignoring demand charge (and comparison) 
% This cell implements the MPC algorithm without explicitly including the
% demand charge.

%
% Zero out the demand charge.
%
  TSDRparams.gd = @(k1,k2,u,w,params) 0;
%
% Call the MPC solver.
%
  [xD, uD, pBarD, costsD] = TSDRMPC(H, w, wnom, TSDRparams);
%
% Redefine the demand charge.
%
  TSDRparams.gd = @(k1,k2,uBig,wBig,params) TSDRgD(k1,k2,uBig,wBig,params);
%
% Compute the cost incurred by the naive policy.
%
  costsD.demandCost = TSDRparams.gd(1,N,uD,w,TSDRparams);
  costsD.totalCost = costsD.totalCost + costsD.demandCost;
%
% Compare the policies with and without including the demand charge. 
%
  TSDRcomparisonPlots(u,costs,uD,costsD,w,5,TSDRparams);
%

%% comparison of MPC and oracle policies 
% This cell compares the power profiles and costs of the MPC and oracle
% policies.

%
% Plot the baseline, MPC and oracle power profiles.
%
  figure(8); clf
  subplot(3,1,1)
  hold on
  stairs(1:48, TSDRparams.pBase, 'r')
  stairs(1:48, u(1,:) + u(3,:) + w(2,:), 'k.-')
  stairs(1:48, uOracle(1,:) + uOracle(3,:) + w(2,:), 'g')
  pNaive = (w(1,:)./TSDRparams.kMain + w(2,:));
  set(gca,'XLim',[0,N],'XTick',0:12:N,'XTickLabel',0:6:round(N/2), ...
    'YLim',[0,50*ceil(max(pNaive)/50)])
  title('Total Power Consumption')
  ylabel('kW')
  legend('Baseline \quad', 'MPC', 'Oracle', 'Location', 'NorthWest')
  set(gca, 'Box', 'on')
%
% Plot the expected and realized cooling loads.
%
  subplot(3,1,2)
  hold on
  stairs(1:48, TSDRparams.wbar(1,:), 'k')
  stairs(1:48, w(1,:), 'r')
  title('Expected and True Cooling Load')
  ylabel('kW$_{th}$')
  set(gca, 'Box', 'on')
  legend('Expected \quad', 'True', 'Location', 'NorthWest')
%
% Plot the expected and realized firm loads.
%
  subplot(3,1,3)
  hold on
  stairs(1:48, TSDRparams.wbar(2,:), 'k')
  stairs(1:48, w(2,:), 'r')
  title('Expected and True Firm Load')
  xlabel('Time (hours)')
  ylabel('kW')
  set(gca, 'Box', 'on')
  legend('Expected \quad', 'True', 'Location', 'NorthWest')
%







