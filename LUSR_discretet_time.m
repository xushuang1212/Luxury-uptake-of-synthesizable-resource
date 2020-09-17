clear all;
close all;

% WITHOUT EVOLUTION WITHOUT MATING
% growth parameters for the leu- strain are labled with 1
% growth parameters for the trp- strain are labled with 2

% leu- strain
Leu = [0.05, 0.1, 0.25, 0.5, 0.75, 1]; Leu = 367.*Leu;
Rmax1 = [0.1459, 0.1912, 0.1989, 0.1859, 0.1649, 0.1708];
Km1 = [51.2176, 63.4227, 33.7455, 32.1082, 26.6437, 29.6791];
Yleu1 = [0.0329, 0.0327, 0.0289, 0.0322, 0.0296, 0.0311]; Yleu1 = 10000000000.*Yleu1;

% trp- strain
Trp = [0.05, 0.1, 0.25, 0.5, 0.75, 1]; Trp = 73.4.*Trp;
Rmax2 = [0.17883, 0.179907, 0.170852, 0.173637, 0.17192, 0.173451];
Km2 = [3.3646, 3.4371, 2.9538, 2.8635, 2.3492, 2.3002];
Ytrp2 = [0.2119, 0.2073, 0.2073, 0.2134, 0.1907, 0.1907]; Ytrp2 = 10000000000.*Ytrp2;

dilution = 2^(-10);
transfer = 3;
transfertime = 300; % time within a transfer (hour)
OD =18500000; % cell number in 1mL when OD = 1
volume = 0.128; % uL, medium volume
Ninitial = volume*OD; % initial population size of both population
freq = [0.1,0.9]; % initial frequency of leu-
death = dilution*transfertime; % death rate per hour

% % Yield -- Unit resource for a species (Yield = (k - N0)/[resource])-number of cells per unit resource
% % Yield should be measured from growth curve data when nutrient is limiting
                  
TRP = [0.05,0.1,0.25,0.5,1];%[0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15];%%%linspace(0.01,1,100);%% % Specific nutrient for population 1 - different concentration of trp (mg/L), fold 0.02 ~ 0.5
LEU = [0.05,0.1,0.25,1];%[0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15];%;%linspace(0.01,1,100);% % Specific nutrient for population 2 - different concentration of leu (mg/L),fold 0.02 ~ 0.25
TRP = 73.4.*TRP; % concetration of nutients (mg/L)
LEU = 367.*LEU;

% LUSR coefficient
LUSR = 0.5; % No luxury uptake
% LUSR = 0.5; % Partial luxury uptake
% LUSR = 1; %Full luxury uptake of synthesizable resource

outputf = zeros(length(LEU),length(TRP),length(freq)); % final frequency of Leu- strain
outputgrid = zeros(length(LEU),length(TRP),length(freq)); % 1 for Leu- strain to increase and -1 to decrease in frequency
   
for i = 1:length(LEU) %leu
      
    for j = 1:length(TRP) %trp
       
        idx1 = find(Trp == TRP(j));
        idx2 = find(Leu == LEU(i));
          
       % yields
        Y1leu = Yleu1(idx1); Y2leu = Y1leu;  
        Y2trp = Ytrp2(idx2); Y1trp = Y2trp; 
         
        for f = 1:length(freq)
             outputscroll = 0;
             
             nutrient = zeros(1,2);
 
             for t = 1:transfer % transfers
                    
                  if t == 1  
                       output = zeros(1,5);
                       Ni = Ninitial*dilution; %binornd(Ninitial,dilution);
                       Ni1 = Ni*freq(f); %binornd(Ni,freq(f)); % initial population of leu-
                       Ni2 = Ni - Ni1; % initial population of trp-        
                  else
                       Ni = output(outputscroll,1)*dilution;
                       Ni1 = Ni*output(outputscroll,4); % initial population of leu-
                       Ni2 = Ni - Ni1; % initial population of trp-
                  end
                  
                  outputscroll = outputscroll + 1;
                  output(outputscroll,1) = Ni; %total population size
                  output(outputscroll,2) = Ni1;% population of Leu-  
                  output(outputscroll,3) = Ni2; %population of trp-
                  output(outputscroll,4) = Ni1/(Ni1+Ni2); % frequency of leu-
                  output(outputscroll,5) = Ni2/(Ni1+Ni2); % frequency of trp-
    
                  
                  leut = LEU(i); nutrient(outputscroll,1) = leut;
                  trpt = TRP(j); nutrient(outputscroll,2) = trpt;
                  
                  for ti = 1:transfertime % within a transfer                      
                      %% growth rate
                         % ynbt = 6.7;leut = leut; trpt = trpt;
                         r1 = Rmax1(idx1)*leut/(Km1(idx1)+leut); r1(r1<=0) = 0;
                         r2 = Rmax2(idx2)*trpt/(Km2(idx2)+trpt); r2(r2<=0) = 0;
                      %% Logistic model of population growth  (ADD CONDITION WHEN STOP GROWING)
                         deltaNi1 = r1*Ni1;
                         Ni1 = Ni1 + deltaNi1;
                         deltaNi2 = r2*Ni2;
                         Ni2 = Ni2 + deltaNi2;
                       %% nutrition update and update growth rate
                       
                        leut = leut - deltaNi1/Y1leu - LUSR * deltaNi2/Y2leu; leut(leut<0) = 0;%Leu
                        trpt = trpt - deltaNi2/Y2trp - LUSR * deltaNi1/Y1trp; trpt(trpt<0) = 0;%Trp
                                          
                       %% output nutrient, population size and frequency over time
                       outputscroll = outputscroll + 1;
                       
                       nutrient(outputscroll,1) = leut; %Leu
                       nutrient(outputscroll,2) = trpt; %Trp                       
                       
                       output(outputscroll,1) = round(Ni1) + round(Ni2); %Total population size
                       output(outputscroll,2) = round(Ni1); %population size of Leu-
                       output(outputscroll,3) = round(Ni2); %population size of Trp-
                       output(outputscroll,4) = output(outputscroll,2)/output(outputscroll,1);%frequency of Leu-
                       output(outputscroll,5) = output(outputscroll,3)/output(outputscroll,1);%frequency of Trp-
                   
                  end
                  %% plot population size
%                   figure(1)
%                   plot(output(:,2),'color', [0 0.5 0]);
%                   hold on
%                   plot(output(:,3),'b');
%                   legend('Leu-','Trp-','Location','best');
%                   xlabel('Time'); ylabel('Population size');
%                   % plot frequency of leu- and trp-
%                   figure(2)
%                   plot(output(:,4),'r'); 
%                   hold on
%                   plot(output(:,5),'k');
%                   xlabel('Time'); ylabel('Frequency of strains');
%                   ylim([0 1]);
%                   legend('Leu-','Trp-','Location','best');
%                   % plot nutrient
%                   figure(3)
%                   plot(nutrient(:,1),'b');
%                   hold on
%                   plot(nutrient(:,2),'r');
%                   legend('Leucine','Tryptophan','Location','best');
%                   xlabel('Time'); ylabel('Nutrient concentration');

                  %% Increase or decrease in frequency for leu-
                  if output(outputscroll,4) > output(1,4)
                       outputgrid(i,j,f) = 1; % leu- increases in frequency
                  elseif output(outputscroll,4) < output(1,4)
                       outputgrid(i,j,f) = -1; % leu- decreases in frequency
                  end
                  
             end
             
%               figure(1)
%               plot(output(:,2),'color', [0 0.5 0]);
%               hold on
%               plot(output(:,3),'b');
%               legend('Leu-','Trp-','Location','best');
%               xlabel('Time'); ylabel('Population size');
%              
%               % plot nutrient
%               figure(3)
%               plot(nutrient(:,1), 'k');
%               hold on
%               plot(nutrient(:,2), 'r');
%               legend('Leucine','Tryptophan','Location','best');
%               xlabel('Time'); ylabel('Nutrient concentration');
        end
    end
end


% outputgrid1 = outputgrid(:,:,1);
% outputgrid2 = outputgrid(:,:,2);
% xy = size(outputgrid1);
% h = zeros(xy(1),xy(2));
% for x = 1:xy(1)
%     for y = 1:xy(2)
%         if outputgrid1(x,y)>0 && outputgrid2(x,y)>0
%             h(x,y) = 3; % leu- wins
%         elseif outputgrid1(x,y)>0 && outputgrid2(x,y)<0
%             h(x,y) = 1; % coexistence
%         elseif outputgrid1(x,y)<0 && outputgrid2(x,y)<0
%             h(x,y) = 2; % trp- wins
%         else
%             h(x,y) = 4;
%         end
%     end
% end
% % h = xlsread('predictionofcoexistence.xlsx',2,'H30:R40');
% figure(4)
% imagesc(h)

