%% ---- DETERMINISTIC, DAY-DEPENDENT TSDR PARAMETERS ---- 
%
%               MPC of Thermal Storage for Demand Response
%               Kevin Kircher, Cornell MAE
%               June 3 '14
%
% This script calculates the physical and economic parameters that depend 
% on the day being simulated. Parameters that may be passed between 
% functions are saved in the struct TSDRparams.
%
% Storage matrices are organized such that the jth column holds the
% quantity at stage j. Columns index time.
%
% SI units are used throughout. Units of kW and kWh refer to electric
% power, while units of kWc and kWhc refer to cooling load. The two are
% related by the coefficient of performance of whatever system is
% providing the cooling:
%
%       COP = (heat removed from space) / (work done by cooling system)
%
% Time is in hours. Money is in dollars.


%% temperature import 
% This cell imports the temperatures at each stage.

%
% Find the hourly temperatures for the simulation day. If the time step
% is half an hour, reshape them into a 1 x 48 array. Use linear
% interpolation.
%
  hourlyTemperatures = cell2mat(...
    hourlyData{iYear,4}(hourlyStartIndex:hourlyEndIndex))';
  if N == 24
    TSDRparams.temperatures = hourlyTemperatures;
  end
  if N == 48
    TSDRparams.temperatures = zeros(1,48);
    TSDRparams.temperatures(1:2:end) = hourlyTemperatures;
    TSDRparams.temperatures(2:2:end) = hourlyTemperatures + ...
      [diff(hourlyTemperatures)/2, 0];
  end
%


%% cooling load mean and variance 
% This cell defines the parameters of the distribution from which the 
% building's cooling load is sampled.

%
% Define the building's cooling load. The model is that at the beginning of
% each stage, a random variable d(k) is realized. This is the new cooling
% load that arrives at stage k. Probabilistically, we treat d(k) as a
% Gaussian random variable with mean dbar(k) and standard deviation
% sigmad(k). 
%
  w1bar = [  161.3605  133.6855  130.8475  123.2115  120.7315  113.8375  111.3920  105.2140  103.1165   97.8575   98.1895  102.4275  106.9825  117.0990  120.3810  134.3820  138.0675  157.8990  162.8290  276.4400  282.9405  303.1560  308.4175  322.9930  327.2290  339.9805  345.8570  334.9540  340.1040  373.5165 375.4255  381.9665  383.1380  383.4355  381.6095  374.7620  370.8830  277.0200  268.0215  225.7720  217.2365  200.6795  195.2230  179.2385  174.4890 162.3725  158.7225  164.6110];
  %sigmaw1 = 1.0e+02*[   1.202123418749999   0.902859200000000   0.848442593750000   0.778247643750000   0.725696206250000   0.649666906250000   0.585155481250000  0.525251393750000   0.533805331250000   0.535118562500000   0.795050981250000   1.658218875000000   2.035236531250000   2.487030393750002  2.105386643750000   1.881159668749999   1.410011656249999   1.403934331250000   1.530926299999999   1.120074187499999   1.218944637500000 1.744910893750001   1.916070656250003   1.722822418749999   1.716615362500000   1.327915231249997   1.229978387500002   2.402477581250000   2.153528893750002   0.568443206249998   0.550453106250000   0.978489331250001   0.730552731250001   0.636141512500001   0.691409387500003   1.243548668749999   1.371321043750003   1.498400812500003   1.591341518750002   1.393689043749999   1.211513143750000   1.249347325000000   1.190296606249998  1.155289456249999   1.099232768750001   1.068404312500000   1.067483406250001   1.056938174999999 ];
  sigmaw1 = 0.1*w1bar;
%
% If the time step is 1 hour, aggregate the half-hour load data. Assume the
% loads are independent from stage to stage, i.e. assume that 
% var(d(k) + d(k-1)) = var(d(k)) + var(d(k-1)) <==> cov(d(k),d(k-1)) = 0.
%
  if N == 24
    w1bar = sum(reshape(w1bar,2,24));
    sigmaw1 = sum(reshape(sigmaw1,2,24));
  end
%


%% firm load mean and variance 
% This cell defines the mean and variance of the fixed load.

%
% Define the building's non-cooling load (lighting, plug loads, etc) means.
% Use the data from Santiago's TRNSYS model, given in kW.
%
  w2bar = zeros(1,N);
  w2bar(1:14) = 40;
  w2bar(15:25) = [47.65 56.41 70.63 85.63 96.25 106.25 113.52 120.4 126.17 131.8 134.7];
  w2bar(26:36) = 137.5;
  w2bar(37:end) = [136.41 135.16 130.63 125.63 117.34 108.6 98.75 88.75 74.38 59.38 49.84 43];
%
% Define the standard deviations.
%
  sigmaw2 = 0.1*w2bar;
%
% If the time step is 1 hour, aggregate the half-hour load data. Assume the
% loads are independent from stage to stage, i.e. assume that 
% var(d(k) + d(k-1)) = var(d(k)) + var(d(k-1)) <==> cov(d(k),d(k-1)) = 0.
%
  if N == 24
    w2bar = sum(reshape(w2bar,2,24));
    sigmaw2 = sum(reshape(sigmaw2,2,24));
  end
%


%% disturbance means and covariance matrices 
% This cell builds the means wbar(k) and covariances Q(k) of the
% disturbance vectors.

%
% Build the means.
%
  TSDRparams.wbar = [w1bar ; w2bar];
%
% Define the correlation coefficient between the cooling load and firm
% load.
%
  rho12 = 0.5*ones(1,N);
%
% Build the covariance matrices and store them in a 3D array.
%
  TSDRparams.Q = zeros(2,2,N);
  TSDRparams.Q(1,1,:) = sigmaw1.^2;
  TSDRparams.Q(2,2,:) = sigmaw2.^2;
  TSDRparams.Q(1,2,:) = (sigmaw1.*sigmaw2).*rho12;
  TSDRparams.Q(2,1,:) = (sigmaw1.*sigmaw2).*rho12;
%


%% ice chiller COPs 
% This section defines the the ice chiller coefficients of performance.

%
% Define the ice chiller coefficients of performance kIce (dimensionless). 
% Note that the relationship between COP and temperature is an output of 
% the TRNSYS building model.
%
  TSDRparams.kIce = min(3.78, 3.7632 - 0.0106*TSDRparams.temperatures);
%


%% main chiller COPs 
% This section defines the main chiller COPs.

%
% Define the main chiller coefficients of performance kMain
% (dimensionless). Note that the relationship between COP and temperature 
% is an output of the TRNSYS building model.
%
  TSDRparams.kMain = min(4.86, 4.837 - 0.0133*TSDRparams.temperatures);
%


%% energy prices 
% This section defines the energy prices.

%
% ---- use the following for a canned Time of Use rate ----
% %
% % Define the time-of-use prices and hours.
% %
%   offPeakPrice = 0.06369;
%   midPeakPrice = 0.19082;
%   peakPrice = 0.23492;
% %
% % Define the energy prices, ce ($/kWh).
% %
%   TSDRparams.ce = offPeakPrice*ones(1,N);
%   if N == 48
%     midPeakStart = 17;
%     peakStart = 25;
%     peakEnd = 36;
%     midPeakEnd = 46;
%     TSDRparams.ce(midPeakStart : midPeakEnd) = midPeakPrice;
%     TSDRparams.ce(peakStart : peakEnd) = peakPrice;
%   end
%   if N == 24
%     midPeakStart = 9;
%     peakStart = 13;
%     peakEnd = 18;
%     midPeakEnd = 23;
%     TSDRparams.ce(midPeakStart : midPeakEnd) = midPeakPrice;
%     TSDRparams.ce(peakStart : peakEnd) = peakPrice;
%   end
%
% ---- use the following for hourly prices based on NYISO LBMPs ----
%
% Define ConEd's Service Classification No. 9 Rate III: General - Large -
% Voluntary Time-of-Day prices. Begin with the flat component, then import
% and add in the NYISO day-ahead LBMPs.
%
  hourlyLBMPs = (hourlyData{iYear,9} ...
    (hourlyStartIndex:hourlyEndIndex)/1000)';
  if N == 24
    TSDRparams.ce = 0.0346*ones(1,N) + hourlyLBMPs;
  end
  if N == 48
    TSDRparams.ce = 0.0346*ones(1,N) + ...
      reshape(repmat(hourlyLBMPs,2,1),1,N);
  end
%


%% demand prices 
% This cell defines the demand prices, which indirectly reflect the
% structure of a monthly demand charge.

%
% Define the incremental $/kW prices for off-peak, mid-peak and peak
% demand, according to ConEd's mandatory hourly pricing rate plan. Use
% these monthly prices in the daily operations in order to encourage load
% shifting to low-price times.
%
  TSDRparams.cdOffPeak = 16.27;
  TSDRparams.cdShoulder = 15.17;
  TSDRparams.cdPeak = 8.10;
% 
% Define the demand prices, cd ($/kW) for each regime (off peak, shoulder
% and peak). Note that the ConEd demand charges are additive: if the 
% monthly peak occurs off-peak, it costs cdOffPeak; if it occurs during a
% mid-peak time (6PM-10PM) it costs cdOffPeak + cdShoulder; and if it occurs 
% during a peak time (8AM-6PM) it costs cdOffPeak + cdShoulder + cdPeak.
% Also note that weekends fall entirely under the off-peak regime.
% 
  TSDRparams.cd = TSDRparams.cdOffPeak*ones(1,N);
  if isbusday(startHour)
    if N == 24
      TSDRparams.shoulderStart = 9;
      TSDRparams.shoulderEnd = 22;
      TSDRparams.cd(TSDRparams.shoulderStart:TSDRparams.shoulderEnd) = ...
        TSDRparams.cdShoulder;
      TSDRparams.peakStart = 9;
      TSDRparams.peakEnd = 18;
      TSDRparams.cd(TSDRparams.peakStart:TSDRparams.peakEnd) = ...
        TSDRparams.cdPeak;
    end
  %  
    if N == 48
      TSDRparams.shoulderStart = 17;
      TSDRparams.shoulderEnd = 44;
      TSDRparams.cd(TSDRparams.shoulderStart:TSDRparams.shoulderEnd) = ...
        TSDRparams.cdShoulder;
      TSDRparams.peakStart = 17;
      TSDRparams.peakEnd = 36;
      TSDRparams.cd(TSDRparams.peakStart:TSDRparams.peakEnd) = ...
        TSDRparams.cdPeak;
    end
  end
%


%% prices of unmet cooling load 
% This cell defines the building occupancy patterns and thermal mass, and
% uses them to compute prices for unmet (or over-met) cooling load.

%
% Define occupancy ratios at each time. Use different schedules for
% business and non-business days. The non-business day occupancies are the
% averages of the Saturday and Sunday/Holiday schedules.
%
  occupancyRatios = ones(1,N);
  if isbusday(startHour)
  %
  % Define the business day occupancy schedule.
  %
    if N == 48
      occupancyRatios(1:12) = 0;
      occupancyRatios(13:14) = 0.1;
      occupancyRatios(15:16) = 0.2;
      occupancyRatios(17:24) = 0.95;
      occupancyRatios(25:26) = 0.5;
      occupancyRatios(27:34) = 0.95;
      occupancyRatios(35:36) = 0.3;
      occupancyRatios(37:44) = 0.1;
      occupancyRatios(45:48) = 0.05;
    end
    if N == 24
      occupancyRatios(1:6) = 0;
      occupancyRatios(7) = 0.1;
      occupancyRatios(8) = 0.2;
      occupancyRatios(9:12) = 0.95;
      occupancyRatios(13) = 0.5;
      occupancyRatios(14:17) = 0.95;
      occupancyRatios(18) = 0.3;
      occupancyRatios(19:22) = 0.1;
      occupancyRatios(23:24) = 0.05;
    end
  else
  %
  % Define the non-business day occupancy schedule.
  %
    if N == 48
      occupancyRatios(1:12) = 0;
      occupancyRatios(13:16) = 0.075;
      occupancyRatios(17:24) = 0.175;
      occupancyRatios(25:34) = 0.075;
      occupancyRatios(35:36) = 0.05;
      occupancyRatios(37:38) = 0.025;
      occupancyRatios(39:48) = 0;
    end
    if N == 24
      occupancyRatios(1:6) = 0;
      occupancyRatios(7:8) = 0.075;
      occupancyRatios(9:12) = 0.175;
      occupancyRatios(13:17) = 0.075;
      occupancyRatios(18) = 0.05;
      occupancyRatios(19) = 0.025;
      occupancyRatios(20:24) = 0;
    end
  end
%
% Calculate the expected occupancy based on the floor area of the building
% and the square feet per person spec. The building has four occupied
% floors, three with area 27ksf and one with area 39ksf (the basement).
%
  floorArea = 120000;
  areaPerPerson = 200;
  maxOccupancy = floorArea/areaPerPerson;
%
% Define the building's thermal mass, Cth (kWhth/K), such that it requires
% Cth kWhth of energy to raise the building's temperature by 1 K.
%
  Cth = 1000*ones(1,N);
%
% Calculate the price of unmet cooling load at each stage, cu 
% ($/(kWhth)^2), such that one kWhth of unmet cooling load incurs a cost of
% cu dollars.
%
% The cost cu is proportional to occupancy (higher cost if more people are
% arround to annoy). It's inversely proportional to the building's thermal
% mass (lower cost if it takes more energy to heat up the building).
%
  TSDRparams.cu = disutilityConstant*maxOccupancy*occupancyRatios./Cth;
%


%% stage and terminal cost functions 
% This section defines anonymous functions to calculate the stage and
% terminal costs. Each component cost function takes a time index and a set
% of column vectors. When these functions are called, the column vectors
% will correspond to the state, control, and/or disturbance vectors at a
% particular stage.

%
% Define the energy cost.
%
  TSDRparams.ge = @(k,uk,wk) TSDRparams.ce(k)*TSDRparams.dt ...
    *(uk(1) + uk(3) + wk(2));
%
% Define the cost of unmet cooling load.
%
% ---- CAUTION the alpha*pos(xk(2)) nerfs conditioning ----
%
%   TSDRparams.gu = @(k,xk) TSDRparams.cu(k)*( xk(2)^2  ...
%     + TSDRparams.alpha*max(0,xk(2)));
%
% Define the stable (but symmetric) cost of unmet load.
%
  TSDRparams.gu = @(k,xk) TSDRparams.cu(k)*( xk(2)^2 );
%
% Define the net stage cost.
%
  TSDRparams.gk = @(k,xk,uk,wk) TSDRparams.ge(k,uk,wk) ...
                 + TSDRparams.gu(k,xk);
%
% Define the tank depletion penalty.
%
  TSDRparams.gTank = @(xN) TSDRparams.ct*pos(TSDRparams.x0(1) - xN(1));
%
% Define the three-tiered demand cost, which is a function of all the
% controls.
%
  TSDRparams.gd = @(k1,k2,uBig,wBig,params) TSDRgD(k1,k2,uBig,wBig,params);
%
% Define the terminal cost.
%
  TSDRparams.gN = @(xN,uBig,wBig) TSDRparams.gTank(xN);
%


%% dynamics function 
% This section defines the dynamics, x(k+1) = f(k,x(k),u(k),w(k)).

%
% Define the dimensions of the problem.
%
  nx = 5;
  nu = 3;
  nw = 2;
%
% Build the nx x nx dynamics matrices, A(k).
%
  TSDRparams.A = zeros(nx,nx,N);
  TSDRparams.A(2,2,:) = 1;
  TSDRparams.A(1,1,:) = TSDRparams.beta;
%
% Build the nx x nu control matrices, B(k).
%
  TSDRparams.B = zeros(nx,nu,N);
  TSDRparams.B(1,1,:) = TSDRparams.dt*TSDRparams.kIce;
  TSDRparams.B(1,2,:) = -TSDRparams.dt./TSDRparams.eta;
  TSDRparams.B(2,2,:) = -TSDRparams.dt;
  TSDRparams.B(2,3,:) = -TSDRparams.kMain*TSDRparams.dt;
  for iTimes = 1 : N
    TSDRparams.B(3:end,:,iTimes) = eye(3);
  end
%
% Build the nx x nw disturbance matrix, G.
%
  TSDRparams.G = zeros(nx,nw);
  TSDRparams.G(2,1) = TSDRparams.dt;
%
% Define the dynamics map, f : R x X x U x R^nw --> X.
%
  TSDRparams.fk = @(k,xk,uk,wk) TSDRparams.A(:,:,k)*xk ...
            + TSDRparams.B(:,:,k)*uk + TSDRparams.G*wk;
%


