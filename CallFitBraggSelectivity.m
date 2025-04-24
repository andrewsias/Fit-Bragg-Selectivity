% Call FitBraggSelectivity

%==========================================================================
% Original disorted holograms.  04 is poor looking data, 08 is hard to fit
%==========================================================================
 load('distortedBraggCurves.mat');

% E05 - strong distortion
% Clean up data 
% [~,maxind] = max(E05(:,1));             % Shift angle to be centered on peak
% dTheta = E05(:,2)-E05(maxind,2);
% 
% minEff = min(E05(:,1));  
% Eff = E05(:,1) - minEff;

% [inE05a,outE05a] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.003,'EffOffset0',0.01,'dPhi20',0.53*2*pi,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi4',1,'FitAlpha',0,'PlotFit',1,'PlotResults',1,'PlotTitle','E05');
% % % 
% [inE05,outE05] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',27.2*10^-6,'n10',0.0068,'EffOffset0',0.03,'dPhi20',0.56*2*pi,'dPhi40',0.05,'FitL',1,'Fitn1',1,...
%     'FitEffOffset',0,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi3',1,'FitdPhi4',1,'FitAlpha',1,'PlotFit',1,'PlotResults',1,'PlotTitle','E05');

% % E06 - strong distortion
% Clean up data 
% [~,maxind] = max(E06(:,1));             % Shift angle to be centered on peak
% dTheta = E06(:,2)-E06(maxind,2);
% 
% minEff = min(E06(:,1));  
% Eff = E06(:,1) - minEff;
% % 
% [inE06,outE06] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',23*10^-6,'n10',0.008,'EffOffset0',0.03,'dPhi20',.49*2*pi,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi4',1,'FitAlpha',0,'PlotFit',1,'PlotResults',1,'PlotTitle','E06');

% E07 - weak distortion
% Clean up data 
% [~,maxind] = max(E07(:,1));             % Shift angle to be centered on peak
% dTheta = E07(:,2)-E07(maxind,2);
% 
% minEff = min(E07(:,1));  
% Eff = E07(:,1) - minEff;
% 
% [inE07,outE07] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',20*10^-6,'n10',0.011,'EffOffset0',0.03,'dPhi20',.01*2*pi,'dPhi40',.07*2*pi,'z04OverL0',0.34,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi3',1,'FitdPhi4',1,'Fitn1',1,'PlotFit',1,'PlotResults',1,'PlotTitle','E07');

% E08 - strong distortion
% Clean up data 
% [~,maxind] = max(E08(:,1));             % Shift angle to be centered on peak
% dTheta = E08(:,2)-E08(maxind,2);
% 
% minEff = min(E08(:,1));  
% Eff = E08(:,1) - minEff;
% 
% % Doesn't fit well, try both Born and not
% 
% [inE08a,outE08a] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',20*10^-6,'n10',0.008,'EffOffset0',0.03,'dPhi20',1.5*pi,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'PlotFit',1,'PlotResults',1,'PlotTitle','E08');
% 
% [inE08b,outE08b] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',outE08a.L,'n10',outE08a.n1,'EffOffset0',outE08a.EffOffset,'dPhi20',outE08a.dPhi2,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi4',1,'PlotFit',1,'PlotResults',1,'PlotTitle','E08');

%==========================================================================
% Exp. 22092202, pre and post cure
% The MATLAB file 22092202_Processed_Data.mat holds the following: preCureData:​
% Item 1 is ‘theta’ which is 8 columns long. Each column corresponds with exposure 01, 02, 03, etc. as detailed in the previous slide.​
% Item 2 is  DE’, which corresponds with theta.​ postCureData:​
% Same as above, but after curing. All columns line up with exposures correctly.
%==========================================================================
 load('22092202_Processed_Data.mat');
% 
% Exp 1 postcure
% dTheta      = postCureData.theta(:,1);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,1);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff(150:210)); 
% maxind      = maxind + 149;
% dTheta      = dTheta - dTheta(maxind);  % Center 

% Polynomial 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,'dPhi20',.3,...
%    'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi3',1,'FitdPhi4',1,'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 1');

%  Fourier
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,'Phi10',.3,...
%    'FitEffOffset',1,'FitThetaOffset',0,'FitPhi1',1,'FitPhi3',1,'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%    'PlotFit',1,'PlotResults',1,'PlotKSpace',0,'PlotTitle','Post cure 1');


% Exp 1 precure
dTheta      = preCureData.theta(:,1);
dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
Eff         = preCureData.DE(:,1);
Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
Eff         = Eff - min(Eff);

[~,maxind]  = max(Eff(150:210)); 
maxind      = maxind + 149;
dTheta      = dTheta - dTheta(maxind);  % Center 

% Polynomial
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,'dPhi20',.9*pi,'dPhi40',-.9*pi,...
%     'FitEffOffset',0,'FitThetaOffset',0,
%     'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 1') ;


% Fourier
[in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
    'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
   'FitEffOffset',0,'FitThetaOffset',0, ...
    'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
     'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
     'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 1') ;

% Exp 2 postcure
% dTheta      = postCureData.theta(:,2);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,2);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff); 
% maxind      = maxind;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% % Polynomial
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,'dPhi20',.3,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi3',1,'FitdPhi4',1,'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 2');
%
% % Fourier
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,'Phi10',.3,...
%     'FitEffOffset',1,'FitThetaOffset',0,'FitPhi1',1,'FitPhi3',1,'FitPhi5',1,...
%     'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%     'PlotFit',1,'PlotResults',1,'PlotKSpace',0,'PlotTitle','Post cure 2');


% % Exp 2 precure
% dTheta      = preCureData.theta(:,2);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = preCureData.DE(:,2);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff(200:220)); 
% maxind      = maxind + 199;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
%
% % Polynomial
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,'dPhi20',.9*pi,'dPhi40',-.9*pi,...
%     'FitEffOffset',0,'FitThetaOffset',0,'FitdPhi2',1,'FitdPhi3',1,'FitdPhi4',1,'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 2');
%
% % Fourier
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
%     'FitEffOffset',0,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 2');
% 
% 
% % Exp 3 postcure
% dTheta      = postCureData.theta(:,3);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,3);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff(150:210)); 
% maxind      = maxind + 149;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% % 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1, ... %     'FitEffOffset',0,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi3',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 3');
% 
% 
% % % Exp 3 precure
% % dTheta      = preCureData.theta(:,3);
% % dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% % Eff         = preCureData.DE(:,3);
% % Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% % Eff         = Eff - min(Eff);
% % 
% % [~,maxind]  = max(Eff(160:200)); 
% % maxind      = maxind + 159;
% % dTheta      = dTheta - dTheta(maxind);  % Center 
% % 
% % [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
% %     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
% %     'FitEffOffset',0,'FitThetaOffset',0,...
% %      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
% %      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
% %      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 3');
% 
% % Exp 4 postcure
% dTheta      = postCureData.theta(:,4);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,4);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff); 
% maxind      = maxind;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,...
%     'FitEffOffset',1,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 4') ;
% 
% % % Exp 4 precure
% dTheta      = preCureData.theta(:,4);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = preCureData.DE(:,4);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff(160:200)); 
% maxind      = maxind + 159;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
%     'FitEffOffset',0,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 4') ;
% 
% 
% % Exp 5 postcure
% dTheta      = postCureData.theta(:,5);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,5);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff); 
% maxind      = maxind;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
%  [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,...
%     'FitEffOffset',1,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 5');
% 
% % % % Exp 5 precure
% dTheta      = preCureData.theta(:,5);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = preCureData.DE(:,5);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff); 
% maxind      = maxind;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.004,'EffOffset0',0.0,...
%     'FitEffOffset',0,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',0,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 5');
% 
% % Exp 6 postcure
% dTheta      = postCureData.theta(:,6);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,6);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff); 
% maxind      = maxind;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,...
%     'FitEffOffset',1,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 6');
% 
% % % Exp 6 precure
% % dTheta      = preCureData.theta(:,6);
% % dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% % Eff         = preCureData.DE(:,6);
% % Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% % Eff         = Eff - min(Eff);
% % 
% % [~,maxind]  = max(Eff(150:175)); 
% % maxind      = maxind+149;
% % dTheta      = dTheta - dTheta(maxind);  % Center 
% % 
% % [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
% %     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
% %     'FitEffOffset',0,'FitThetaOffset',0,...
% %      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
% %      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
% %      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 6');
% 
% % Exp 7 postcure
% dTheta      = postCureData.theta(:,7);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,7);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff); 
% maxind      = maxind;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,...
%     'FitEffOffset',0,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 7');
% % 
% % % % Exp 7 precure
% dTheta      = preCureData.theta(:,7);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = preCureData.DE(:,7);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = max(Eff(190:210)); 
% maxind      = maxind+189;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',1,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
%     'FitEffOffset',1,'FitThetaOffset',0,...
%      'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 7');
% 
% % Exp 8 postcure
% dTheta      = postCureData.theta(:,8);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = postCureData.DE(:,8);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = min(Eff(160:190));                     % Note MIN
% maxind      = maxind+159;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',30*10^-6,'n10',0.006,'EffOffset0',0.1,...
%     'FitEffOffset',0,'FitThetaOffset',0,...
%     'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Post cure 8');
% 
% % % % Exp 8 precure
% dTheta      = preCureData.theta(:,8);
% dTheta      = dTheta(~isnan(dTheta));       % Remove NaNs at end
% Eff         = preCureData.DE(:,8);
% Eff         = Eff(~isnan(Eff));             % Remove NaNs at end
% Eff         = Eff - min(Eff);
% 
% [~,maxind]  = min(Eff(170:190));        % note MIN
% maxind      = maxind+169;
% dTheta      = dTheta - dTheta(maxind);  % Center 
% 
% [in1,out1] = FitBraggSelectivity(dTheta,Eff,'Born',0,'RefWrtTheta',-14.6625,'ObjWrtTheta',14.6625, ...
%     'L0',25*10^-6,'n10',0.008,'EffOffset0',0.0,...
%     'FitEffOffset',0,'FitThetaOffset',0, ...
%     'FitPhi1',1,'FitPhi2',1,'FitPhi3',1,'FitPhi4',1,'FitPhi5',1,...
%      'FitA1',1,'FitA2',1,'FitA3',1,'FitA4',1,'FitA5',0,'FitA6',0,...
%      'PlotFit',1,'PlotResults',1,'PlotTitle','Pre cure 8');
% 
% 
% %-------
% % Summary of results pre cure and post cure fitting with and without amplitude variations
% %--------
% LPreAmp     = [23.38 25.84 24.08 28.97 25.51 24.02 28.03 27.32];
% LPostAmp    = [22.4  26.25 22.76 25.6  28.25 24.13 28.86 29.86];
% nAvgPreAmp  = [7.41   4.78  7.15  7.7  7.67   8.37  8.42  4.76];
% nAvgPostAmp = [8.47   6.89  7.84 8.34  7.62   6.68  4.89  5.26];
% 
% LPreNoAmp   = [30.86 34.52];
% LPostNoAmp  = [23.65 28.63];
% n1PreNoAmp  = [6.23  5.99];
% n1PostNoAmp = [8.58  6.38];
% 
% figure;
% hold off;
% plot(1:8,LPreAmp,'b+','MarkerSize',12);hold on
% plot(1:8,LPostAmp,'b.','MarkerSize',12);
% plot(1:2,LPreNoAmp,'r+','MarkerSize',12);
% plot(1:2,LPostNoAmp,'r.','MarkerSize',12);
% xlim([0,9]);
% ylim([0,35]);
% set(gcf,'color','w');                           % Plot background white
% set(gca,'FontSize',18);
% set(gca,'FontName','Arial');
% xlabel('Experiment #');
% title('Change in thickness upon curing');
% ylabel('L [{\mu}m]');
% legend('Pre (Amp & phase)','Post (Amp & phase)','Pre (Phase only)','Post (Phase only)','Location','SE');
% 
% figure;
% hold off;
% plot(1:8,nAvgPreAmp,'b+','MarkerSize',12);hold on
% plot(1:8,nAvgPostAmp,'b.','MarkerSize',12);
% plot(1:2,n1PreNoAmp,'r+','MarkerSize',12);
% plot(1:2,n1PostNoAmp,'r.','MarkerSize',12);
% xlim([0,9]);
% ylim([0,10]);
% set(gcf,'color','w');                           % Plot background white
% set(gca,'FontSize',18);
% set(gca,'FontName','Arial');
% xlabel('Experiment #');
% title('Change in average n_1 upon curing');
% ylabel('<n_1>');
% legend('Pre (Amp & phase)','Post (Amp & phase)','Pre (Phase only)','Post (Phase only)','Location','SE');
% 
