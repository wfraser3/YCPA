winstyle = 'docked';
%winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15; %Simulation time
f = 500e12; %Frequency
lambda = c_c/f; %Wavelength

xMax{1} = 20e-6; %Setting up region
nx{1} = 200; %Steps in the x direction
ny{1} = 0.75*nx{1}; %Steps in the y direction


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0; %Permeability

epi{1} = ones(nx{1},ny{1})*c_eps_0; %Permittivity
inclusions = input('Type 0 for grating or 1 for cool structure: ');
if(inclusions == 0)
    epi{1}(125:150,55:95)= c_eps_0*11.3;
    epi{1}(125:150,1:41) = c_eps_0*11.3;
    epi{1}(125:150,110:150) = c_eps_0*11.3;
elseif(inclusions == 1)
    epi{1}(50:60,20:25) = c_eps_0*11.3;
    epi{1}(65:75,20:25) = c_eps_0*11.3;
    epi{1}(80:90,20:25) = c_eps_0*11.3;
    epi{1}(95:105,20:25) = c_eps_0*11.3;
    epi{1}(110:120,20:25) = c_eps_0*11.3;
    epi{1}(125:135,20:25) = c_eps_0*11.3;
    epi{1}(140:150,20:25) = c_eps_0*11.3;
    
    epi{1}(50:60,125:130) = c_eps_0*11.3;
    epi{1}(65:75,125:130) = c_eps_0*11.3;
    epi{1}(80:90,125:130) = c_eps_0*11.3;
    epi{1}(95:105,125:130) = c_eps_0*11.3;
    epi{1}(110:120,125:130) = c_eps_0*11.3;
    epi{1}(125:135,125:130) = c_eps_0*11.3;
    epi{1}(140:150,125:130) = c_eps_0*11.3;
    
    epi{1}(140:150,30:35) = c_eps_0*11.3;
    epi{1}(140:150,40:45) = c_eps_0*11.3;
    epi{1}(140:150,50:55) = c_eps_0*11.3;
    epi{1}(140:150,60:65) = c_eps_0*11.3;
    epi{1}(140:150,70:75) = c_eps_0*11.3;
    epi{1}(140:150,80:85) = c_eps_0*11.3;
    epi{1}(140:150,90:95) = c_eps_0*11.3;
    epi{1}(140:150,100:105) = c_eps_0*11.3;
    epi{1}(140:150,110:115) = c_eps_0*11.3;
end

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});

dx = xMax{1}/nx{1}; %Step size
dt = 0.25*dx/c_c; %Time step
nSteps = round(tSim/dt*2); %Simulation length
yMax = ny{1}*dx; %Y dimension size
nsteps_lamda = lambda/dx; %Number of lambdas

movie = 1; %Lets make a movie
Plot.off = 0; %Plot toggle
%A bunch of plotting stuff:
Plot.pl = 0; 
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

%Boundary conditions
bc{1}.NumS = 1;
bc{1}.s(1).xpos = nx{1}/(4) + 1;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
%Wave parameters
mag = 1;
phi = 0;
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
%st = 15e-15;
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; %Boundary conditions

bc{1}.s(2).xpos = nx{1}/(4) + 1;
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;
% mag = -1/c_eta_0;
%Wave parameters
mag = 5;
phi = pi;
omega = f*2*pi;
betap = 0;
t0 = 45e-15;
%st = 15e-15;
st = -0.05;
s = 0;
y0 = yMax/2;
sty = 2*lambda;
bc{1}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'}; %Boundary conditions

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg






