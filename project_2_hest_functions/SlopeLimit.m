function psi = SlopeLimit(a,b,type)
% function psi = SlopeLimit(a,b,type,c,M,h)
% Define slope limiter function based on type. c is used in one limiter
% M,h is used for TVB limiting
N = length(a); zero = zeros(N,1);

% No slope 
if (type==0) psi = zeros(N,1); end
% minmod limiter
if (type==1) psi = minmod([a b]); end
% MUSCL limiter
if (type==2) psi = minmod([(a+b)/2 2*a 2*b]); end
return