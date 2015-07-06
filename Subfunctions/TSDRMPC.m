function [x, u, pBar, costs] = TSDRMPC(H, w, wnom, params)
%TSDRMPC computes the MPC policy for the TSDR problem.
%
% This function implements a receding horizon, single-stage roll-out
% algorithm with nominal values substituted in for random variables in the
% TES4SR problem. This is known as model-predictive control (MPC).
%
% Starting from k = 0, the MPC algorithm is
%
% 1. Update the current nominal disturbance by setting wbar(k) = w(k)
%
% 2. Solve the kth deterministic optimal control subproblem. Begin by  
%    defining the effective horizon, L = min(H,N-k), where H is the MPC 
%    horizon.
% 
%      piMPC = arg min_pi  JMPC(k,x(k),wbar)
%              subject to  u(k+j) in U
%                          x(k+j+1) in X
%                          x(k+j+1) = f(k+j,x(k+j),u(k+j),wbar(k+j))
% 
%    where the constraints apply for all j in {0, ..., L-1}. This results
%    in the MPC policy
% 
%     piMPC = [uMPC(k), uMPC(k+1), ..., uMPC(k+L-1)]
% 
%    Note the use of the nominal disturbances in the cost function
%    and dynamics.
% 
%    The cost function JMPC in the optimal control subproblem is
% 
%      JMPC(k,x0,wbar) = gN(k+L) 
%                        + sum_{j=k}^{k+L-1} gk(j,x(j),u(j),wbar(j))
%       
% 3. Implement only the first control, u(k) = uMPC(k), and let the state
%    evolve
%
% 4. If k = N-1, stop; otherwise, increment k and go to 2
%
%
% ---- Inputs ----
%
%   H                   The MPC horizon, a scalar in [1,N].
%
%   w                   = [w(0), w(1), ..., w(N-1)] is the nw x nTimes 
%                       matrix of disturbances, the jth column of which 
%                       holds the disturbance vector at stage j.
%
%   wnom                A nw x nTimes matrix, the jth column of which holds
%                       the nominal disturbance vector at stage j.
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
%   pBar                = [pBar(0), pBar(1), ..., pBar(N-1)] is the 3 x 
%                       nTimes matrix of running peak demands (one for
%                       each tier of the demand charge).
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
% Initialize the state vectors x(0), x(1), ..., x(N) and the control 
% vectors u(0), u(1), ..., u(N-1).
%
  x = zeros(nx,N+1);
  x(:,1) = x0;
  u = zeros(nu,N);
%
% Declare cost storage vectors. These will keep track of the following, at
% every stage: the cost-to-go predicted by the MPC algorithm, the energy 
% cost, the curtailment cost, the cost of unmet cooling load, and the cost
% of noncompliance with deployment calls.
%
  predictedCostsToGo = zeros(1,N);
  energyCosts = zeros(1,N);
  DRcosts = zeros(1,N);
  disutilityCosts = zeros(1,N);
%
% ---- MPC loop ----
%
% Note that the state vectors on paper are indexed from 0 to N or from 
% 0 to N-1 for everything else. Since Matlab can't index from zero, this
% loop uses the shifted time index kp = k+1, so that the state vectors
% are defined for kp = 1, 2, ..., N+1. Everything else is defined for
% kp = 1, 2, ..., N.
%
  for kp = 1 : N
  %
  % ---- Step 2 ----
  %
  % Calculate the horizon for the subproblem.
  %
    L = min(H,N-kp+1);
  %
  % Enter CVX scope to define and solve the optimal control subproblem.
  %
    cvx_begin quiet
    %
    % Declare the decision variable, piMPC = [uMPC(kp), ..., uMPC(kp+L)].
    % Also declare the intermediate decision variable 
    % xMPC = [x(k), ..., x(k+L)] to be resolved by the dynamics 
    % constraint.
    %
      variable piMPC(nu,L);
      variable xMPC(nx,L+1);
    %
    % Define the cost function.
    %
      if kp == 1
        JMPC = gN(xMPC(:,L+1)) + gd(kp,kp+L-1,piMPC,wnom,params);
      else
        JMPC = gN(xMPC(:,L+1)) + gd(kp,kp+L-1,piMPC,wnom(:,kp:kp+L-1),params);
      end
      for j = 1 : L
        JMPC = JMPC + gk(kp+j-1,xMPC(:,j),piMPC(:,j),wnom(:,kp+j-1));
      end
    %
    % Tell CVX to minimize the cost function.
    %
      minimize JMPC
    %
    % Tell CVX that the constraints follow.
    %
      subject to
    %
    % Each control vector must be in the (time-invariant) feasible control
    % space.
    %
      repmat(uMin,1,L) <= piMPC <= repmat(uMax,1,L);
      for j = 1 : L
        abs(piMPC(:,j) - xMPC(3:end,j)) <= params.duMax;
      end
    %
    % The states, which are functions of the controls and the nominal
    % disturbances, must obey the dynamics at all stages.
    %
      xMPC(:,1) == x(:,kp);
      for j = 1 : L
        xMPC(:,j+1) == fk(kp+j-1,xMPC(:,j),piMPC(:,j),wnom(:,kp+j-1));
      end
    %
    % The first state must be in the time-invariant feasible state space.
    % Note that the second state is unconstrained.
    %
      x1min <= xMPC(1,:) <= x1max;
    %
    % Solve the problem and exit CVX scope.
    %
    cvx_end
  %
  % Store the first MPC control vector and the predicted cost-to-go under
  % the full MPC policy.
  %
    u(:,kp) = piMPC(:,1);
    predictedCostsToGo(kp) = cvx_optval;
  %
  % Calculate and store the stage costs.
  %
    energyCosts(kp) = ge(kp,u(:,kp),w(:,kp));
    DRcosts(kp) = gDR(kp,u(:,kp),w(:,kp));
    disutilityCosts(kp) = gu(kp,x(:,kp));
  %
  % Print an update note to the workspace.
  %
    fprintf('The current MPC iteration is %d.\n', kp);
  %
  % ---- Step 5 ----
  %
  % Propagate the state using the MPC control law and the true disturbance.
  %
    x(:,kp+1) = fk(kp,x(:,kp),u(:,kp),w(:,kp));
  %
  % ---- Step 6 ----
  %
  % Update the running peak demand over all hours.
  %
    params.pBar(1,kp+1) = max(params.pBar(1,kp), [1 0 1]*u(:,kp) + w(2,kp));
  %
  % Update the running peak demand for the shoulder tier.
  %
    if kp >= params.shoulderStart && kp <= params.shoulderEnd
      params.pBar(2,kp+1) = max(params.pBar(2,kp), [1 0 1]*u(:,kp) + w(2,kp));
    end
    if kp <= params.shoulderStart - 1 || kp >= params.shoulderEnd+1
      params.pBar(2,kp+1) = params.pBar(2,kp);
    end
  %
  % Update the running peak demand for the peak tier.
  %
    if kp >= params.peakStart && kp <= params.peakEnd
      params.pBar(3,kp+1) = max(params.pBar(3,kp), [1 0 1]*u(:,kp) + w(2,kp));
    end
    if kp <= params.peakStart - 1 || kp >= params.peakEnd+1
      params.pBar(3,kp+1) = params.pBar(3,kp);
    end  
  %
  % End of simulation loop.
  %
  end
%
% Prepare the running peak power matrix for export.
%
  pBar = params.pBar;
%
% Calculate the demand cost, terminal cost and total cost.
%
  demandCost = gd(1,N,u,w,params);
  tankCost = gN(x(:,kp+1));
  totalCost = sum(energyCosts + DRcosts + disutilityCosts) + ...
    demandCost + tankCost;
%
% Build the output struct.
%
  costs.predictedCostsToGo = predictedCostsToGo;
  costs.energyCosts = energyCosts;
  costs.DRcosts = DRcosts;
  costs.disutilityCosts = disutilityCosts;
  costs.demandCost = demandCost;
  costs.tankCost = tankCost;
  costs.totalCost = totalCost;
%
% End of function TESMPC.
%
  end

