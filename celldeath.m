% Nathan Liu, John Cocjin, Gabrien Clark, Navin Pathak
% BIOE 446: Computational Modeling Lab
% Dr. Amina Qutub

% Cell Death Model
% Assumptions:
% - no cell proliferation
% - media is well-mixed above plated cells in wells
% - no direct cell-cell (gap junction) interactions
% - no chemical flow in/out of cells while healthy
tic
clc
clear all;
close all;
statemat = [];
casp3mat = [];


% 1 Well of a well plate - 10 k cells
n = 9;
n_x = sqrt(n);
n_y = sqrt(n);
concmat = [];

% time scale - milliseconds
dt = 8640;
tend_h = 24;
tend_m = tend_h*24;
tend_s = tend_m*60;
tv = 1:dt:tend_s;

timemat = [];

% Radius of hepatocyte: if V = 3.4*10^-9 cm^3, r = 20.8153 um
r = 20.8153*10^-6;
d = 2*r;
rad = num2cell(r*ones(1,n));

% position data
xmap = (r:d:r+(n_x-1)*d);
ymap = (r:d:r+(n_y-1)*d);

xp = zeros(n_x,n_y);
yp = zeros(n_x,n_y);

for i = 1:1:n_x
    for j = 1:1:n_y
        xp(i,j) = xmap(i);
        yp(i,j) = ymap(j);
    end
end

xp1 = reshape(xp,1,n_x*n_y);
yp1 = reshape(yp,1,n_x*n_y);

xpos = num2cell(xp1);
ypos = num2cell(yp1);

% State of cell (1 = live, 0 = apoptotic)
state = num2cell(ones(1,n));

% TNFalpha concentration
% base conc.
tnf = 20; % using relative 20:1 tnfalpha:tnfr ratio - actual: 25 ng/mL
tnfa = zeros(n_x,n_y);


% Concentration Gradient
xmin = 1;
xmax = n_x;
xcenter = (xmin+xmax)/2;

ymin = 1;
ymax = n_y;
ycenter = (ymin+ymax)/2;

rmax = sqrt((xmax-xcenter)^2+(ymax-ycenter)^2);

for i = 1:1:n_x
    for j = 1:1:n_y
        % distance from center of cells
        dist_center = sqrt((i-xcenter)^2 + (j - ycenter)^2);
        % min radius = 0 conc = 20
        % max radius = n_x conc = 0
        
        tnfa(i,j) = tnf*((dist_center - rmax)/rmax)^2;
    end
end

tnfalpha = num2cell(reshape(tnfa,1,n_x*n_y));

complex1 = num2cell(ones(1,n));
% half-life TNF - 18 mins
halflife_m = 18;
halflife_s = halflife_m*60;


% TNFR1 receptor concentration
tnfr1 = num2cell(ones(1,n));

% Caspase 8/10 concentration
casp8 = num2cell(zeros(1,n));

% Caspase 3 concentration
casp3 = num2cell(zeros(1,n));

casp8i = num2cell(ones(1,n));

casp3i = num2cell(10*ones(1,n));


% Initialization of structure used to hold values
cells = struct('xpos',xpos,'ypos',ypos,'state',state,'radius',rad,...
    'tnfr1',tnfr1,'complex1',complex1,'tnfalpha',tnfalpha,...
    'casp8',casp8,'casp3',casp3,'casp8i',casp8i,'casp3i',casp3i,...
    'conc',concmat,'time',timemat);

% For loop for iteration at each time step
for it = 2:length(tv)
    for ic = 1:length(cells)
        tspan = (tv(it-1):tv(it));
        
        % calculate change in state/ROS/TNFalpha, run func.
        [y,t] = cdeath(cells(ic),tspan);
        if find(y(:,5)>=0.8*cell2mat(casp3i(ic)),1) ~= 0
            cells(ic).state = 0;
        end
        
        cells(ic).conc = [cells(ic).conc;y];
        cells(ic).time = [cells(ic).time;t];
        
        cells(ic).tnfalpha = y(end,1);
        cells(ic).tnfr1 = y(end,2);
        cells(ic).complex1 = y(end,3);
        cells(ic).casp8 = y(end,4);
        cells(ic).casp3 = y(end,5);
        cells(ic).casp8i = y(end,6);
        cells(ic).casp3i = y(end,7);
    end
    
   statemat = [statemat;[cells.state]];
   casp3mat = [casp3mat;[cells.casp3]];
end
% Gabrien wants position,radius,caspase3,state
toc