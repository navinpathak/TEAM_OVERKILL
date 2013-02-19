% inputs:
%  cells - structure containing cells information
%  tspan - time span
% outputs:
%  y - matrix of concentrations as a function of time
%  t - vector of time

function [y,t] = cdeath(cells,tspan)
k = [.001 .05 .08 .1];       % rates of reactions as defined in project
                            % k(1) <- binding of TNF-alpha and TNFR1
                            % k(2) <- unbinding of TNF-alpha - TNFR1 complex
                            % k(3) <- activation rate of caspase 8/10
                            % k(4) <- activation rate of caspase 3

%k(1) = .05;
%k(3) = .8;

y0 = [cells.tnfalpha cells.tnfr1 cells.complex1 cells.casp8 cells.casp3 ...
    cells.casp8i cells.casp3i];

[t,y] = ode45(@ckin,tspan,y0,[],k);

end

function dx = ckin(t,y,k)

% compound y(1) <- TNF-alpha
% compound y(2) <- TNFR1
% compound y(3) <- TNF-alpha - TNFR1
% compound y(4) <- Caspase 8
% compound y(5) <- Caspase 3
% compound y(6) <- Inactive caspase 8
% compound y(7) <- Inactive caspase 3

t12 = 1080;

dx(1) = -k(1)*y(1)*y(2) + k(2)*y(3) - log(2)/t12*y(1);
dx(2) = -k(1)*y(1)*y(2) + k(2)*y(3);
dx(3) = +k(1)*y(1)*y(2) - k(2)*y(3);
dx(4) = k(3)*y(3)*y(6);
dx(5) = k(4)*y(4)*y(7);
dx(6) = -k(3)*y(3)*y(6);
dx(7) = -k(4)*y(4)*y(7);

dx = dx';

end