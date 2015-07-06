function [x, u, costs] = TSDROLOC(w, params)
%TSDROLOC computes the OLOC policy for the TSDR problem.
%
% This function implements an open-loop optimal control algorithm in the
% TSDR problem. This function solves the single deterministic optimal 
% control problem:
% 
%    piOLOC = arg min_pi  JOLOC(x(0),...,x(N),w)
%             subject to  u(k) in U(x(k))
%                         x(k+1) in X
%                         x(k+1) = f(k,x(k),u(k),w(k))
%
% where the constraints apply for all k in {0, ..., N-1}. This results in
% the OLOC policy
% 
%    piOLOC = [uOLOC(0), uOLOC(1), ..., uOLOC(N-1)]
% 
% The cost function in the optimal control problem is
% 
%    JOLOC(k,x0,w(0),...,w(N-1)) = gN(x(N),u(0),...,u(N-1),w(0),...,w(N-1)) ...
%                                + sum_{k=0}^{N-1} gk(k,x(k),u(k),w(k))
%
% If the argument w passed to this function is the matrix of true disturbance 
% vectors, then this function calculates the optimal policy under the 
% oracle information pattern (prescient knowledge of all disturbances at
% stage 0).
%
% If the arguments w is not the matrix of true disturbance vectors, then
% this function will return a suboptimal policy.
%
%
% ---- Inputs ----
%
%   w                   = [w(0), w(1), ..., w(N-1)] is the nw x nTimes 
%                       matrix of disturbances, the jth column of which 
%                       holds the disturbance vector at stage j.
%
% params is a struct containing:
%
%   gk                  An anonymous function handle for the stage cost.
%                       Inputs to gk are the stage k, the state x(k), the
%                       control u(k), and the disturbance w(k). The output
%                       from gk is the scalar stage cost.
%
%   gN                  An anonymous function handle for the terminal cost.
%                       The input to gN is the terminal state x(N). The
%                       output from gN is the scalar terminal cost.
%
%   gd                  An anonymous function handle for the demand cost.
%                       The inputs to gd are the nu x N matrix of controls
%                       and the nw x N matrix of disturbances.
%
%   f                   An anonymous function handle for the state
%                       dynamics. The inputs to f are the stage k, the
%                       state x(k), the control u(k) and the disturbance
%                       w(k). The output of f is the state at stage k+1,
%                       x(k+1).
%
%   uMin                The nu x 1 lower bound on the feasible controls.
%
%   uMax                The nu x 1 upper bound on the feasible controls.
%
%   x1Min               The scalar lower bound on x1.
%
%   x1Max               The scalar upper bound on x1.
%
%
% ---- Outputs ----
%
%   x                   = [x(0), x(1), ..., x(N)] is the nx x (nTimes+1)
%                       matrix of states, the jth column of which stores 
%                       the state at stage j.
%
%   u                   = [u(0), u(1), ..., u(N-1)] is the nu x nTimes
%                       matrix of controls, the jth column of which stores
%                       the control applied at stage j.
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
%   demandCost          The scalar peak demand cost.
%
%   tankCost            The scalar cost of tank depletion.
%
%   totalCost           The scalar total cost.

%
% Import the relevant data from the input struct.
%
  x0 = params.x0;
  uMin = params.uMin;
  uMax = params.uMax;
  x1min = params.x1min;
  x1max = params.x1max;
%
% Import the cost and dynamics functions.
%
  fk = params.fk;
  ge = params.ge;
  gDR = params.gDR;
  gd = params.gd;
  gu = params.gu;
  gk = params.gk;
  gN = params.gN;
%
% Import the problem dimensions.
%
  [nx,nu,N] = size(params.B);
%
% ---- OLOC solution ----
%
% Note that the state vectors on paper are indexed from 0 to N and the
% controls and disturbances are indexed from 0 to N-1. Since Matlab can't 
% index from zero, this loop uses the shifted time index kp = k+1, so that
% the state vectors are defined for kp = 1, 2, ..., N+1. Everything else is
% defined for kp = 1, 2, ..., N.
%
  cvx_begin quiet
  %
  % Declare the decision variable, piOLOC = [uOLOC(0), ..., uOLOC(N-1)].
  % Also declare the intermediate decision variable 
  % xOLOC = [x(0), ..., x(N)] to be resolved by the dynamics 
  % constraint.
  %
    variable piOLOC(nu,N);
    variable xOLOC(nx,N+1);
  %
  % Define the cost function.
  %
    JOLOC = gN(xOLOC(:,N+1),piOLOC,w) + gd(1,N,piOLOC,w,params);
    for kp = 1 : N
      JOLOC = JOLOC + gk(kp,xOLOC(:,kp),piOLOC(:,kp),w(:,kp));
    end
  %
  % Tell CVX to minimize the cost function.
  %
    minimize JOLOC
  %
  % Tell CVX that the constraints follow.
  %
    subject to
  %
  % Each control vector must be in the (time-invariant) feasible control
  % space.
  %
    repmat(uMin,1,N) <= piOLOC <= repmat(uMax,1,N);
    for kp = 1 : N
      abs(piOLOC(:,kp) - xOLOC(3:end,kp)) <= params.duMax;
    end
  %
  % The states, which are functions of the controls and the nominal
  % disturbances, must obey the dynamics at all stages.
  %
    xOLOC(:,1) == x0;
    for kp = 1 : N
      xOLOC(:,kp+1) == fk(kp,xOLOC(:,kp),piOLOC(:,kp),w(:,kp));
    end
  %
  % The first element of each state vector must be in the time-invariant 
  % feasible state space. Note that the second element of the state is 
  % unconstrained.
  %
    x1min <= xOLOC(1,:) <= x1max;
  %
  % Solve the problem and exit CVX scope.
  %
  cvx_end
%
% Store the OLOC control and state vectors.
%
  u = piOLOC;
  x = xOLOC;
%
% Calculate the stage costs.
%
  energyCosts = zeros(1,N);
  DRcosts = zeros(1,N);
  disutilityCosts = zeros(1,N);
  for kp = 1 : N
    energyCosts(kp) = ge(kp,u(:,kp),w(:,kp));
    DRcosts(kp) = gDR(kp,u(:,kp),w(:,kp));
    disutilityCosts(kp) = gu(kp,x(:,kp));
  end
%
% Calculate the demand cost, terminal cost and total cost.
%
  demandCost = gd(1,N,u,w,params);
  tankCost = gN(x(:,kp+1),u,w);
  totalCost = sum(energyCosts + DRcosts + disutilityCosts) ...
    + demandCost + tankCost;
%
% Build the output struct.
%
  costs.energyCosts = energyCosts;
  costs.DRcosts = DRcosts;
  costs.disutilityCosts = disutilityCosts;
  costs.demandCost = demandCost;
  costs.tankCost = tankCost;
  costs.totalCost = totalCost;
%
% End of function TSDROLOC.
%
  end

