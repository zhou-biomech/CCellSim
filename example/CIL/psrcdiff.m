function varargout = psrcdiff(varargin) % point source diff
% Postma M, Van Haastert P J M. A diffusion–translocation model for 
% gradient sensing by chemotactic cells[J]. Biophysical journal, 2001, 
% 81(3): 1314-1323.
%
%            cb      rb      cosh( (rb-r)/lambda )         inf
%  c(r,t) = ---- ---------- ----------------------- - cb * sum    tau_n cos(n*pi*r/rb) * exp(-t/tau)
%            kd    lambda      sinh( rb/lambda )           n=-inf
%
% varargin = r           varargout = cinf
% varargin = [r, t]      varargout = [c, cinf]
% varargin = [C, R, f]   varargout = r0

D = 10;          % diffusion rate  % 100
kd = 0.05;        % degradation rate

rb = 200;         % boundary condition
cb = 0.00041145;

n = reshape(-100:100, 1, 1, 201);
lambda = sqrt(D/kd);
gamma = cb/kd * rb/lambda / sinh(rb/lambda);
tau = 1./(kd + n.^2 * pi^2 *D / rb^2);

switch nargin
    case 1
        r = varargin{:};
        cinf =  gamma * cosh((rb-r)/lambda);
      % grad = -gamma * sinh((rb-r)/lambda)/lambda;
        varargout{1} = cinf;
    case 2
        [r, t] = varargin{:};
        cinf =  gamma * cosh((rb-r)/lambda);
        c = cinf - cb * sum(tau .* cos(n * pi .* r/rb) .* exp(-t/tau), 3);
        c = max(c,0);
        varargout{1} = c;
        varargout{2} = cinf;
    case 3
        [C, R, f] = varargin{:};
        r0 = lambda*log((exp(rb/lambda)*( ...
             C^2*f^2*exp((2*R)/lambda) ...
             -(2+2*f+f^2)*gamma^2*exp((2*R)/lambda) ...
             +(1+f)*gamma^2*exp((4*R)/lambda) ... 
             +(1+f)*gamma^2)^(1/2) ...
             -C*f*exp((R + rb)/lambda))/(gamma*(exp((2*R)/lambda) ...
             + f*exp((2*R)/lambda) - 1)));
        varargout{1} = r0;
end


% -------------------------------------------------------------------------
% syms gamma rb r lambda R C f x
% c = @(y) gamma * cosh((rb-y)/lambda);
% solve(c(x-R) - c(x+R) - f * (c(x+R)+C), x)
