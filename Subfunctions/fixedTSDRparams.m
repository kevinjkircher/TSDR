%% ---- FIXED TSDR PARAMETERS ---- 
%
%               MPC of Thermal Storage for Demand Response
%               Kevin Kircher, Cornell MAE
%               June 3, 2014
%
% This script generates the fixed physical and economic parameters for 
% the TSDR building simulation. Parameters that may be passed between
% functions are saved in the struct CVXTSDRparams.
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


%% timing 
% This section defines the control horizon and step size.

%
% Define the control horizon, T (hours).
%
  T = 24;
%
% Define the step size, dt (hours), and the number of time steps, N.
%
  TSDRparams.dt = 0.5;
  N = T / TSDRparams.dt;
  TSDRparams.N = N;
%


%% ice chiller parameters 
% This section defines the physical parameters of the ice chiller.

%
% Define the ice chiller capacity u1max (kW), the minimum ice chiller 
% power u1min (kW), the ice tank storage capacity x1max (kWhc), and the
% maximum ice extraction per stage (kWhc). Also define the maximum ramp
% rates du1max (kW) and du2max (kWhc).
%
  rampFraction = 0.75;
  kWperTon = 3.51685;
  u1max = 94;
  u1min = 0*u1max;
  du1max = rampFraction*u1max;
  TSDRparams.x1min = 0;
  TSDRparams.x1max = 500*kWperTon; % NOTE: 800 ton-hr is Santiago's capacity
  u2min = 0;
  u2max = 145*kWperTon;
  du2max = 1*u2max;
%
% Define the ice dissipation rates beta (dimensionless) and the ice 
% extraction efficiencies eta (dimensionless).
%
  TSDRparams.beta = 0.98*ones(1,N);
  TSDRparams.eta = 0.9*ones(1,N);
%


%% main chiller parameters 
% This section defines the physical parameters of the main chiller.

%
% Define the main chiller capacity u3max (kW) and the minimum ice chiller
% power u3min (kW). Also define the maximum ramp rate, du3max (kW).
%
  u3max = 73;
  u3min = 0*u3max;
  du3max = rampFraction*u3max;
%
% For convenience, collect all the control constraints into vectors.
%
  TSDRparams.uMax = [u1max ; u2max ; u3max];
  TSDRparams.uMin = [u1min ; u2min ; u3min];
  TSDRparams.duMax = [du1max ; du2max ; du3max];
%

