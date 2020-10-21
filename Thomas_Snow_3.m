function ModVar = Thomas_Snow_3(Data, ParVal, InState, Prog)
% function ModVar = Thomas_SSim(Data, ParVal, InState, Prog)
% 10/05/05 Hoshin V. Gupta (Original for HyMod01)
% 01/31/06 Guillermo Martinez (Modified for Thomas Model)
% 01/08/08 Guillermo Martinez (Paper Modification)
% INPUTS
%   Data = Data Structure
%   Parval = Structure containing Parameters for Simulation
%   InState = Initial State structure
%   Prog = Structure of program variables
% OUTPUTS
%   ModVar = Model Computed variables structure - final values
%=========================================================================
%--(0)--Preallocate space for variables
Period = Prog.SimPer;       % Simulation Period
Nmonth = length(Period);
W = zeros(1,Nmonth);
Y = zeros(1,Nmonth);    
S = zeros(1,Nmonth);
G = zeros(1,Nmonth);
Q = zeros(1,Nmonth);
E = zeros(1,Nmonth);
GR = zeros(1,Nmonth);
MTMIN = zeros(1,Nmonth);
MTMAX = zeros(1,Nmonth);

MPP = zeros(1,Nmonth);    
MPE = zeros(1,Nmonth);

%--Snow components
%SP = zeros(1,Nmonth);
SN = zeros(1,Nmonth);     
SM = zeros(1,Nmonth);       

%--(1)--Initialize variables - flag = 0 for Det & flag = 1 for Stoch
flag = 0;
if Prog.SimOpt == 2; flag = 1; end;

%Data
%PP  = Data.PP;                % Precipitation during the month
%PET = Data.PET;               % PET during the month

%Snow modelling data
TP    = Data.PP;              % Total precipitation during the month
PET   = Data.PET;             % PET during the month
TMIN  = Data.tmin;            % Mininum Temperature
TMAX  = Data.tmax;            % Maximum Temperature    
%default_tsnow   = Data.tsnow; % Threshold where only snow is present
default_train   = 3.3;        % Threshold where only rain is present
default_meltmax = 0.5;        % Maximum melting rate
default_deltatsnow = 0.0;     % Delta tsnow


% Parameters stochastic simulation
SigPP = flag*ParVal.SigPP;
SigPE = flag*ParVal.SigPE;    
SigW  = flag*ParVal.SigW;    
%SigY  = flag*ParVal.SigY;
    
% Parameters Thomas model (abcd)   
a = ParVal.a;   % Runoff occurence before the soil is saturated [0<=a<=1] 
b = ParVal.b;   % Limit sum evapotransporation and soil moisture storage
c = ParVal.c;   % Distribution factor between groundwater recharge and direct runoff [0 <= c <= 1]
d = ParVal.d;   % Groundwater residence time

% Initial states
   
%Select method to calculate the initial states
if Prog.calS == 1
   S0 = ParVal.e;  % Calibration SCE
else
   S0 = InState.S; % Default value
end

if Prog.calG == 1
   G0 = ParVal.f;  % Calibration SCE
else
   G0 = InState.G; % Default value
end       

% Parameters for snow modeling
if Prog.caltrain == 1
   train = ParVal.g;  % Calibration SCE
else
   train = default_train; % Default value
end  

% Parameters for snow modeling
if Prog.caldelta == 1
   deltatosnow=ParVal.j;
   tsnow = train-deltatosnow;  % Calibration SCE
else
   tsnow = default_deltatsnow; % Default value
end  

% Parameters for snow modeling
if Prog.calmeltmax == 1
   meltmax = ParVal.i;  % Calibration SCE
else
   meltmax = default_meltmax; % Default value
end  

[PP,SP] = snowpartition(TP,TMIN,TMAX,tsnow,train,Prog.snowopt);  

SN0=0; %Initial Storage of Snow

%--(2)--Run Model Simulation TIME Loop
      
for i = 1:Nmonth;
% .S = Soil moisture storage at the end of the month
% .E = Actual evaporation
% .G = Groundwater storage at the end of the month        
        
   Month = Period(i);
   MPP(i) = PP(Month);
   MTMIN(i) = TMIN(Month);    
   MTMAX(i) = TMAX(Month);
   if PP(Month)>0; 
      PPeps = SigPP*randn*log10(PP(Month)); 
      MPP(i) = max( 10^(log10(PP(Month)) + PPeps), 0 ); 
   end;
   
   %Snow state calculation
   if i > 1
      SN(i) = SN(i-1) + SP(i);
   else
      SN(i) = SN0 + SP(i);     
   end
   
   %SM(i) = max((T(i)-tsnow)/(train-tsnow)*meltmax*SN(i),0);
   SM(i) = meltsnow(SN(i),meltmax,MTMIN(i),MTMAX(i),tsnow,train,Prog.snowopt);

   SN(i) = SN(i) - SM(i);
    
%    if(SN(i) < 10)
%       SM(i) = SM(i)+SN(i);
%       SN(i) = 0;
%    end  
       
   PEeps = SigPE*randn*PET(Month);
   MPE(i)= max(PET(Month)+PEeps,0);               
                
   if i > 1
      W(i) = MPP(i) + S(i-1) + SM(i);
   elseif 1 == i
      W(i) = MPP(i) + S0;                    
   end
        
   Weps = SigW * rand * W(i);
   W(i) = max(W(i)+Weps ,0);              
        
   %Y(i) = ((W(i) - e) + b)/(2*a)-((((W(i) - e) + b)/(2*a))^2 - (W(i) - e)*b/a )^(1/2);
   Y(i) = ((W(i)) + b)/(2*a)-((((W(i)) + b)/(2*a))^2 - (W(i))*b/a )^(1/2);
        
   Yeps = SigW * rand * Y(i);
        
   Y(i) = max(Y(i)+Yeps ,0);        
        
   S(i) = Y(i)*exp(-MPE(i)/b);
        
   E(i) = Y(i) - S(i);
        
   if i >1
      G(i) = (G(i-1)+c*(W(i)-Y(i)))/(1+d);
   elseif 1 == i
      G(i) = (G0+c*(W(i)-Y(i)))/(1+d);
   end

   Q(i) = (1-c)*(W(i)-Y(i))+d*G(i);
        
   GR(i) = c*(W(i)-Y(i));
                
end;
    
%--(3)-- Send model data out
ModVar.W  = W';           % Model computed available water state contents
ModVar.Y  = Y';           % Model computed evapotranspiration opportuniy state contents
ModVar.PP = MPP;          % Model computed precip
ModVar.PE = MPE;          % Model computed pot evap
%Solving error for dates with a different start month
ModVar.TMIN = MTMIN;          % Model computed pot evap
ModVar.TMAX = MTMAX;          % Model computed pot evap
ModVar.S  = S;            % Model computed moisture storage
ModVar.Q  = Q;            % Model computed runoff
ModVar.G  = G;            % Model computed groundwater flow   
ModVar.E  = E;            % Model computed actual evaporation
ModVar.GR = GR;           % Model computed groundwater recharge
%-- Variables from snow component
ModVar.PP = PP(Period);         % Model computed precipitation as precipitation
ModVar.SP = SP(Period);         % Model computed precipitation as snow
ModVar.SM = SM;                 % Model computed snowmelt
ModVar.SN = SN;                 % Model computed snowpack
ModVar.TP = TP(Period);         % Model total precipitation
ModVar.train = train;           % Model total precipitation
ModVar.tsnow = tsnow;           % Model total precipitation
    
%-- End of function Thomas_Snow_3.m