function [RAIN,SNOW] = snowpartition(TOTAL,TMIN,TMAX,tsnow,train,optmodel)
%Guillermo Martinez 09/28/2006
%Revision HWR677 12/03/2007
%Empirical model to separate monthly rain from snow using temperature
%
%INPUTS:
%TOTAL:     Array - Total Precipitation
%TMIN:      Array - Average Minimum Temperature
%TMAX:      Array - Average Maximum Temperature
%tsnow:     Threshold where only snow is present
%train:     Threshold where only rain is present
%optmodel:  Option used to split between rain and snow
%           1 - Only use TMIN
%           2 - Only use TMAX
%           3 - Combination of TMIN and TMAX
%
%OUTPUTS:
%SNOW:      Array - Total Precipitation as Snow
%RAIN:      Array - Total Precipitation as Rain

SNOW=zeros(length(TOTAL),1);
RAIN=zeros(length(TOTAL),1);

%Select only snow, intermediate or only rain for each case
switch optmodel
   case 1;
      [allrain]    = find(TMIN > train);
      [rainorsnow] = find(TMIN <= train & TMIN >= tsnow);
      [allsnow]    = find( TMIN < tsnow);
      
      SNOW(rainorsnow) = TOTAL(rainorsnow).*((train - TMIN(rainorsnow)))/(train - tsnow);
   case 2;
      [allrain]    = find(TMAX > train);
      [rainorsnow] = find(TMAX <= train & TMAX >= tsnow);
      [allsnow]    = find( TMAX < tsnow);
      SNOW(rainorsnow) = TOTAL(rainorsnow).*((train - TMAX(rainorsnow)))/(train - tsnow);
   case 3;  
      allrain=zeros(lengtT(TOTAL),1); rowrain=1;
      allsnow=zeros(lengtT(TOTAL),1); rowsnow=1;
      rainorsnow=zeros(lengtT(TOTAL),1);rowrasn=1;
      
      for i=1:length(TOTAL)
         if (TMIN(i)<tsnow)
            allsnow(rowsnow)=i;                    
            rowsnow = rowsnow +1;            
         elseif (TMAX(i)>train)
            allrain(rowrain)=i;
            rowrain = rowrain +1;            
         else
            rainorsnow(rowrasn)=i;
            rowrasn = rowrasn +1;
         end         
      end
      allsnow = nonzeros(allsnow);
      allrain = nonzeros(allrain);
      rainorsnow = nonzeros(rainorsnow);      
      TAVG=TMAX-TMIN;
      SNOW(rainorsnow) = TOTAL(rainorsnow).*((train - TAVG(rainorsnow)))/(train - tsnow);
      %keyboard
end
  
%Assign values depending according to the ranges of temperature  
if (allrain)
   RAIN(allrain) = TOTAL(allrain);
   SNOW(allrain) = 0;
end

if (rainorsnow)
   RAIN(rainorsnow) = TOTAL(rainorsnow) - SNOW(rainorsnow);
end

if (allsnow)
   RAIN(allsnow) = 0;
   SNOW(allsnow) = TOTAL(allsnow);
end
