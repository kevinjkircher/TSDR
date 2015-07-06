function dCost = TSDRgD(k1,k2,uBig,wBig,params)
%TSDRcD computes the demand charge increase for the TSDR problem.
%
% This function computes the increase in the demand charge between stages
% k1 and k2 >= k1 in the TSDR problem.
%
% ---- Inputs ----
%
%   k1                  the start time index, a scalar in [0,N-1].
%
%   k2                  the end time index, a scalar in [k1,N-1].
%
%   uBig                =[u1(k1),...,u1(k2); ...; u3(k1),...,u3(k2)], the 
%                       3 x k2-k1+1 matrix of control vectors.
%
%   wBig                =[w1(k1),...w1(k2); w2(k1),...,w2(k2)], the 
%                       2 x k2-k1+1 matrix of disturbance vectors.
%
%   params              a parameter struct containing:
%
%       cdOffPeak           the off-peak price ($/kW).
%
%       cdShoulder          the shoulder price ($/kW).
%
%       cdPeak              the peak price ($/kW).
%
%       shoulderStart       the scalar index of the shoulder start time.
%
%       shoulderEnd         the scalar index of the shoulder end time.
%
%       peakStart           the scalar index of the peak start time.
%
%       peakEnd             the scalar index of the peak end time.
%
%       pBar                =[pBar(k1),...,pBar(k2)], the 3 x (k2-k1+1)
%                           matrix of running peak power values in each of 
%                           the three demand tiers.
%
% ---- Outputs ----
%
%   dCost                   the scalar demand cost ($).

%
% Make sure the time arguments are specified correctly. If not, return a
% cost of zero and print an error message.
%
  if k2 < k1
    dCost = 0;
    fprintf('Error: demand cost function called with k2 < k1.')
    return
  end
%

%% tier 1: all hours 
% This cell computes the first tier of the demand charge, which covers 
% all hours.

%
% Calculate the off-peak charge.
%
  dCostOffPeak = params.cdOffPeak*pos(max([1 0 1]*uBig + wBig(2,:)) ...
    - params.pBar(1,k1));
%


%% tier 2: shoulder + peak 
% This cell computes the second tier of the demand charge, which covers 
% both shoulder and peak.

%
% Case 1: k1,k2 < shoulderStart
%
  if k2 < params.shoulderStart
    dCostShoulder = 0;
  end
%
% Case 2: k1 <= shoulderStart <= k2 <= shoulderEnd
%
  if k1 <= params.shoulderStart && params.shoulderStart <= k2 ...
      && k2 <= params.shoulderEnd
    kStart = 1 + params.shoulderStart - k1;
    dCostShoulder = params.cdShoulder*pos(max([1 0 1]*uBig(:,kStart:end) ...
      + wBig(2,kStart:end)) - params.pBar(2,k1));
  end
%
% Case 3: k1 <= shoulderStart < shoulderEnd < k2
%
  if k1 <= params.shoulderStart && params.shoulderEnd < k2
    kStart = 1 + params.shoulderStart - k1;
    kEnd = 1 + params.shoulderEnd - k1;
    dCostShoulder = params.cdShoulder*pos(max([1 0 1]*uBig(:,kStart:kEnd) ...
      + wBig(2,kStart:kEnd)) - params.pBar(2,k1));
  end
%
% Case 4: shoulderStart < k1 <= k2 <= shoulderEnd
%
  if params.shoulderStart < k1 && k2 <= params.shoulderEnd
    dCostShoulder = params.cdShoulder*pos(max([1 0 1]*uBig + wBig) ...
      - params.pBar(2,k1));
  end
%
% Case 5: shoulderStart < k1 <= shoulderEnd < k2
%
  if params.shoulderStart < k1 && k1 <= params.shoulderEnd ...
      && params.shoulderEnd < k2
    kEnd = params.shoulderEnd - k1 + 1;
    dCostShoulder = params.cdShoulder*pos(max([1 0 1]*uBig(:,1:kEnd) ...
      + wBig(2,1:kEnd)) - params.pBar(2,k1));
  end    
%
% Case 6: shoulderStart < shoulderEnd < k1 <= k2
%
  if params.shoulderEnd < k1
    dCostShoulder = 0;
  end
%


%% tier 3: peak 
% This cell computes the third tier of the demand charge, which covers the
% peak.

%
% Case 1: k1,k2 < peakStart
%
  if k2 < params.peakStart
    dCostPeak = 0;
  end
%
% Case 2: k1 <= peakStart <= k2 <= peakEnd
%
  if k1 <= params.peakStart && params.peakStart <= k2 ...
      && k2 <= params.peakEnd
    kStart = 1 + params.peakStart - k1;
    dCostPeak = params.cdPeak*pos(max([1 0 1]*uBig(:,kStart:end) ...
      + wBig(2,kStart:end)) - params.pBar(3,k1));
  end
%
% Case 3: k1 <= peakStart < peakEnd < k2
%
  if k1 <= params.peakStart && params.peakEnd < k2
    kStart = 1 + params.peakStart - k1;
    kEnd = 1 + params.peakEnd - k1;
    dCostPeak = params.cdPeak*pos(max([1 0 1]*uBig(:,kStart:kEnd) ...
      + wBig(2,kStart:kEnd)) - params.pBar(3,k1));
  end
%
% Case 4: peakStart < k1 <= k2 <= peakEnd
%
  if params.peakStart < k1 && k2 <= params.peakEnd
    dCostPeak = params.cdPeak*pos(max([1 0 1]*uBig + wBig) - params.pBar(3,k1));
  end
%
% Case 5: peakStart < k1 <= peakEnd < k2
%
  if params.peakStart < k1 && k1 <= params.peakEnd ...
      && params.peakEnd < k2
    kEnd = params.peakEnd - k1 + 1;
    dCostPeak = params.cdPeak*pos(max([1 0 1]*uBig(:,1:kEnd) ...
      + wBig(2,1:kEnd)) - params.pBar(3,k1));
  end    
%
% Case 6: peakStart < peakEnd < k1 <= k2
%
  if params.peakEnd < k1
    dCostPeak = 0;
  end
%


%% total demand cost 
% This cell adds up the components of the demand cost.

%
% Add the off-peak, shoulder and peak contributions.
%
  dCost = dCostOffPeak + dCostShoulder + dCostPeak;
%

%
% End of function TSDRgD.
%
  end

