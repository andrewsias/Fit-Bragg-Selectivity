%==========================================================================
% Fit Bragg Selectivity
%
% Input parameters: 
%   DeltaThetaExt   External angle of read beam relative to Bragg matched angle [deg]
%   EffVsDeltaTheta Measured diffraction efficiency at each angle
%   Optional: Name,value pairs as defined by dictionary in parser
%
% Outputs as fields of the "out" struct. Upper and lower confidence bounds also reported.
%   n1              Peak-to-mean index variation at read wavelength.  
%   L               Material thickness [m]
%   alpha           Exponential grating decay (>0) of grating amplitude in z [1/m]
%   dPhi[n]         [n]th polynomial phase distortion [rad]
%   z04             Location of symmetric nulls in quartric distortion [m] 0<=z04<=L
%   Phi[n]          [n]th Fourier series phase distortion [rad]
%
% History
% RRM         Oct 1 2022    Originate
%==========================================================================
function [in,out,n1,L,GoF] = FitBraggSelectivity(DeltaThetaExt,EffVsDeltaTheta,varargin)

[in,out]  = ParseFBSInputs(varargin);

out = FitExpToTheory(in,out,DeltaThetaExt,EffVsDeltaTheta);

if in.PlotFit
    PlotResults(DeltaThetaExt,EffVsDeltaTheta,in,out);
    if in.PlotResults, PlotFitResults(in,out); end
    if in.PlotPhi,     PlotPhi(in,out);        end
    if in.PlotA,       PlotA(in,out);          end
    if in.PlotKSpace,  PlotKSpace(in,out);     end
end

% Add special outputs for Andrew (isolate var for parfor):
n1  = out.n1;
L   = out.L;
GoF = out.gof.rsquare;

end % FitBraggSelectivity

%==========================================================================
% Parse inputs
%==========================================================================
function [in,out] = ParseFBSInputs(VariablesPassed)

nano  = 10^-9;                      % units.  All calculations MKS.
%micro = 10^-6;
%milli = 10^-3;

p = inputParser;                    % Parser structure

%---Parser validation functions. 
isPosScalar  = @(x) isnumeric(x) && isscalar(x) && (x >= 0);     % Scalar > 0

%------------------------------Parameter dictionary------------------------
%                       PARAM             DEFAULT  VALIDITY      DEFINITION
%--------------------------------------------------------------------------
%---EXPERIMENTAL PARAMETERS (passed thru to BraggTransmission function)
addParameter(p,          'Born',            false, @isscalar   ); % Assume Born (weak) approximation
addParameter(p,   'PowerNormal',            false, @isscalar   ); % TRUE: power is measured normal to surface with an overfilled 
                                                                  % detector (so Sz).  FALSE: power is measured normal to k or the 
                                                                  % detector area is larger than the beam

addParameter(p,   'RefWrtTheta',            -22.5, @isscalar   ); % Reference write beam angle from surface normal [deg]
addParameter(p,   'ObjWrtTheta',             22.5, @isscalar   ); % Object write beam angle from surface normal [deg]
addParameter(p,    'lambdaWrt0',         405*nano,  isPosScalar); % Vacuum write optical wavelength
addParameter(p,    'lambdaRed0',         635*nano,  isPosScalar); % Vacuum  read optical wavelength
addParameter(p,          'nWrt',              1.5,  isPosScalar); % Refractive index at write wavelength
addParameter(p,          'nRed',              1.5,  isPosScalar); % Refractive index at read wavelength

%---PLOT CONTROLS
addParameter(p,       'PlotFit',             true, @isscalar   ); % Plot results?
addParameter(p,    'PlotKSpace',            false, @isscalar   ); % Show k space insert on fit results?
addParameter(p,       'PlotPhi',             true, @isscalar   ); % Show distortion subplot
addParameter(p,       'PlotA',               true, @isscalar   ); % Show attenuation subplot
addParameter(p,   'PlotResults',             true, @isscalar   ); % Show table of fit results on plot?
addParameter(p,     'PlotTitle',               '', @ischar     ); % Title for plots

%---WHAT TO FIT  
addParameter(p,          'FitL',            true, @isscalar    ); % Vertical shift of data to compensate for background light
addParameter(p,         'Fitn1',            true, @isscalar    ); % Vertical shift of data to compensate for background light
addParameter(p,  'FitEffOffset',           false, @isscalar    ); % Vertical shift of data to compensate for background light
addParameter(p,'FitThetaOffset',           false, @isscalar    ); % Horizontal shift of data to compensate for Bragg mismatch
addParameter(p,      'Fitalpha',           false, @isscalar    ); % Exponential grating decay (>0) of grating amplitude in z [1/m]
%-Polynomial phase distortion
addParameter(p,      'FitdPhi2',           false, @isscalar    ); % Peak quadratic phase distortion [rad]
addParameter(p,      'FitdPhi3',           false, @isscalar    ); % Peak cubic phase distortion [rad]
addParameter(p,      'FitdPhi4',           false, @isscalar    ); % Coefficient of quartic phase distortion [rad]
%-Fourier series phase distortion
addParameter(p,       'FitPhi1',           false, @isscalar    ); % Phase [rad] of of first  harmonic (Period = 2 L / 1)
addParameter(p,       'FitPhi2',           false, @isscalar    ); % Phase [rad] of of second harmonic (Period = 2 L / 2)
addParameter(p,       'FitPhi3',           false, @isscalar    ); % Phase [rad] of of third  harmonic (Period = 2 L / 3)
addParameter(p,       'FitPhi4',           false, @isscalar    ); % Phase [rad] of of fourth harmonic (Period = 2 L / 4)
addParameter(p,       'FitPhi5',           false, @isscalar    ); % Phase [rad] of of fifth  harmonic (Period = 2 L / 5)
addParameter(p,       'FitPhi6',           false, @isscalar    ); % Phase [rad] of of sixth  harmonic (Period = 2 L / 6)
%-Fourier series ampliutde
addParameter(p,         'FitA1',           false, @isscalar    ); % Amlitude [relative] of of first  harmonic (Period = 2 L / 1)
addParameter(p,         'FitA2',           false, @isscalar    ); % Amlitude [relative] of of second harmonic (Period = 2 L / 2)
addParameter(p,         'FitA3',           false, @isscalar    ); % Amlitude [relative] of of third  harmonic (Period = 2 L / 3)
addParameter(p,         'FitA4',           false, @isscalar    ); % Amlitude [relative] of of fourth harmonic (Period = 2 L / 4)
addParameter(p,         'FitA5',           false, @isscalar    ); % Amlitude [relative] of of fifth  harmonic (Period = 2 L / 5)
addParameter(p,         'FitA6',           false, @isscalar    ); % Amlitude [relative] of of sixth  harmonic (Period = 2 L / 6)

%---START POINTS: Optimization is faster from good initial guesses.  Use
% simpler fits to estimate values for more complex (more parameter) fits
addParameter(p,            'L0',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'n10',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,    'EffOffset0',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,  'ThetaOffset0',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,        'alpha0',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,        'dPhi20',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,        'dPhi30',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,        'dPhi40',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,     'z04OverL0',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,         'Phi10',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,         'Phi20',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,         'Phi30',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,         'Phi40',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,         'Phi50',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,         'Phi60',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'A10',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'A20',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'A30',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'A40',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'A50',             NaN, @isscalar    ); % Starting point for fit
addParameter(p,           'A60',             NaN, @isscalar    ); % Starting point for fit

%---FIT COEFFICIENT SCALING: Advanced option.
% Nonlinear least squares fitting (as implemented
% by Matlab) works best when all the coefficients are order unity.  Scaling
% factors are used to pass scaled versions of variables far from unity
% (like thickness), then rescale the results.  Consider modifying these
% when coefficients times the default scaling factor are far from unity.
% Typically this improves both speed and final results of fit.
addParameter(p,          'ScaleL',          10^5, isPosScalar  );
addParameter(p,         'Scalen1',          10^2, isPosScalar  );
addParameter(p,  'ScaleEffOffset',             1, isPosScalar  ); 
addParameter(p,'ScaleThetaOffset',             1, isPosScalar  );
addParameter(p,      'Scalealpha',         10^-5, isPosScalar  );
addParameter(p,      'ScaledPhi2',             1, isPosScalar  );
addParameter(p,      'ScaledPhi3',             1, isPosScalar  );
addParameter(p,      'ScaledPhi4',             1, isPosScalar  );
addParameter(p,   'Scalez04OverL',             1, isPosScalar  );
addParameter(p,       'ScalePhi1',             1, isPosScalar  );
addParameter(p,       'ScalePhi2',             1, isPosScalar  );
addParameter(p,       'ScalePhi3',             1, isPosScalar  );
addParameter(p,       'ScalePhi4',             1, isPosScalar  );
addParameter(p,       'ScalePhi5',             1, isPosScalar  );
addParameter(p,       'ScalePhi6',             1, isPosScalar  );
addParameter(p,         'ScaleA1',             1, isPosScalar  );
addParameter(p,         'ScaleA2',             1, isPosScalar  );
addParameter(p,         'ScaleA3',             1, isPosScalar  );
addParameter(p,         'ScaleA4',             1, isPosScalar  );
addParameter(p,         'ScaleA5',             1, isPosScalar  );
addParameter(p,         'ScaleA6',             1, isPosScalar  );

%--------------------------------------------------------------------------
% Parse
%--------------------------------------------------------------------------
parse(p,VariablesPassed{:});                        % Parse inputs into struct p
in = p.Results;                                     % Short hand notation 

%--------------------------------------------------------------------------
% Random setup and input massaging
%--------------------------------------------------------------------------
% Names of all parameters that need to be passed thru to BraggTransmission
in.PassThruParams = {'Born','PowerNormal','RefWrtTheta','ObjWrtTheta','lambdaWrt0','lambdaRed0','nWrt','nRed'};

% Names of fit parameters given an initial value and not fit.  ALSO pass
% these through to BraggTransmission
AllFitCoefNames     = {'FitL','Fitn1','FitEffOffset','FitThetaOffset','Fitalpha',...
    'FitdPhi2','FitdPhi3','FitdPhi4','FitPhi1','FitPhi2','FitPhi3','FitPhi4','FitPhi5','FitPhi6',...
    'FitA1','FitA2','FitA3','FitA4','FitA5','FitA6'};
in.PassThruCoeffs   = '';

for ic = 1:length(AllFitCoefNames)
    % Not fitting coefficient but initial value specified AND not one of
    % the Offset parameters that are not passed to BraggTransmission
    if ~in.(AllFitCoefNames{ic}) && ~isnan(in.([AllFitCoefNames{ic}(4:end),'0'])) ...
            && ~contains(AllFitCoefNames{ic},'Offset')
        in.PassThruCoeffs = [in.PassThruCoeffs, {AllFitCoefNames{ic}(4:end)}]; % Strip "Fit" to get pass thru parameter name. Add to list
    end
end

% Assign fixed values of L and n1 if not fitting
out.RunDate = datetime('today');         % Initililze out
if ~in.FitL  && ~isnan(in.L0),  out.L  = in.L0;  out.LUpperCI  = in.L0;  end
if ~in.Fitn1 && ~isnan(in.n10), out.n1 = in.n10; out.n1UpperCI = in.n10; end

end % ParseFBSInputs

%==========================================================================
% Fit data
%==========================================================================
function out = FitExpToTheory(in,out,DeltaThetaExt,EffVsDeltaTheta)

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
OptionsExp          = '';           % Initialize expression strings used to construct fittype
EffOffsetExp        = '';
ThetaOffsetExp      = 'x';          % If not fitting offset in angle (x variable of fit)

StartPoint          = [];           % Initilize vectors of options
Lower               = [];
Upper               = [];

%--------------------------------------------------------------------------
% Create string expressions and fit option vectors that will be used in the
% fitting call.  The vectors of options (starting point and bounds) must be
% in the order that matlab chooses to arrange the fitting parameters.  This
% turns out to be upper case then lower case, both alphabetical.  So the
% code below MUST be arranged in this same manner.  Any new variables added
% has to be inserted in the right order or the code will break.
%--------------------------------------------------------------------------
if in.FitA1         % An are first in list
    OptionsExp      = [OptionsExp, ',''A1'',A1/',num2str(in.ScaleA1)];      
    if isnan(in.A10)
        StartPoint  = [StartPoint,  0*in.ScaleA1];
    else
        StartPoint  = [StartPoint,  in.A10*in.ScaleA1];   
    end
    Lower           = [Lower,   -1*in.ScaleA1];
    Upper           = [Upper,    1*in.ScaleA1];                       % Amplitude is relative to 1
end

if in.FitA2        
    OptionsExp      = [OptionsExp, ',''A2'',A2/',num2str(in.ScaleA2)];      
    if isnan(in.A20)
        StartPoint  = [StartPoint,  0*in.ScaleA2];
    else
        StartPoint  = [StartPoint,  in.A20*in.ScaleA2];   
    end
    Lower           = [Lower,   -1*in.ScaleA2];
    Upper           = [Upper,    1*in.ScaleA2];                       % Amplitude is relative to 1
end

if in.FitA3        
    OptionsExp      = [OptionsExp, ',''A3'',A3/',num2str(in.ScaleA3)];      
    if isnan(in.A30)
        StartPoint  = [StartPoint,  0*in.ScaleA3];
    else
        StartPoint  = [StartPoint,  in.A30*in.ScaleA3];   
    end
    Lower           = [Lower,   -1*in.ScaleA3];
    Upper           = [Upper,    1*in.ScaleA3];                       % Amplitude is relative to 1
end

if in.FitA4        
    OptionsExp      = [OptionsExp, ',''A4'',A4/',num2str(in.ScaleA4)];      
    if isnan(in.A40)
        StartPoint  = [StartPoint,  0*in.ScaleA4];
    else
        StartPoint  = [StartPoint,  in.A40*in.ScaleA4];   
    end
    Lower           = [Lower,   -1*in.ScaleA4];
    Upper           = [Upper,    1*in.ScaleA4];                       % Amplitude is relative to 1
end

if in.FitA5        
    OptionsExp      = [OptionsExp, ',''A5'',A5/',num2str(in.ScaleA5)];      
    if isnan(in.A50)
        StartPoint  = [StartPoint,  0*in.ScaleA5];
    else
        StartPoint  = [StartPoint,  in.A50*in.ScaleA5];   
    end
    Lower           = [Lower,   -1*in.ScaleA5];
    Upper           = [Upper,    1*in.ScaleA5];                       % Amplitude is relative to 1
end

if in.FitA6        
    OptionsExp      = [OptionsExp, ',''A6'',A6/',num2str(in.ScaleA6)];      
    if isnan(in.A60)
        StartPoint  = [StartPoint,  0*in.ScaleA6];
    else
        StartPoint  = [StartPoint,  in.A60*in.ScaleA6];   
    end
    Lower           = [Lower,   -1*in.ScaleA6];
    Upper           = [Upper,    1*in.ScaleA6];                       % Amplitude is relative to 1
end

if in.FitEffOffset  % Starts with upper case E
    EffOffsetExp    = ['EffOffset/',num2str(in.ScaleEffOffset),' + '];     % Include efficiency offset
    if isnan(in.EffOffset0)
        StartPoint  = [StartPoint,  0*in.ScaleEffOffset];
    else
        StartPoint  = [StartPoint,  in.EffOffset0*in.ScaleEffOffset];
    end
    Lower           = [Lower,    -0.2*in.ScaleEffOffset];
    Upper           = [Upper,    +0.2*in.ScaleEffOffset];                  % +/- 20% offset
end

if in.FitL          % Starts with upper case L                            
    OptionsExp      = [OptionsExp, ',''L'',L/',num2str(in.ScaleL)];    
    if isnan(in.L0)
        StartPoint  = [StartPoint,  20E-6*in.ScaleL];                      % Start at L = 20 um
    else
        StartPoint  = [StartPoint,  in.L0*in.ScaleL];
    end
    Lower           = [Lower,           0*in.ScaleL];
    Upper           = [Upper,        1E-3*in.ScaleL];                      % Max of L = 1 mm
end

if in.FitPhi1
    OptionsExp      = [OptionsExp, ',''Phi1'',Phi1/',num2str(in.ScalePhi1)];      
    if isnan(in.Phi10)
        StartPoint  = [StartPoint,  0*in.ScalePhi1];
    else
        StartPoint  = [StartPoint,  in.Phi10*in.ScalePhi1];   
    end
    Lower           = [Lower,   -4*pi*in.ScalePhi1];
    Upper           = [Upper,    4*pi*in.ScalePhi1];                       % Two waves of distortion
end

if in.FitPhi2
    OptionsExp      = [OptionsExp, ',''Phi2'',Phi2/',num2str(in.ScalePhi2)];      
    if isnan(in.Phi20)
        StartPoint  = [StartPoint,  0*in.ScalePhi2];
    else
        StartPoint  = [StartPoint,  in.Phi20*in.ScalePhi2];   
    end
    Lower           = [Lower,   -4*pi*in.ScalePhi2];
    Upper           = [Upper,    4*pi*in.ScalePhi2];                       % Two waves of distortion
end

if in.FitPhi3
    OptionsExp      = [OptionsExp, ',''Phi3'',Phi3/',num2str(in.ScalePhi3)];      
    if isnan(in.Phi30)
        StartPoint  = [StartPoint,  0*in.ScalePhi3];
    else
        StartPoint  = [StartPoint,  in.Phi30*in.ScalePhi3];   
    end
    Lower           = [Lower,   -4*pi*in.ScalePhi3];
    Upper           = [Upper,    4*pi*in.ScalePhi3];                       % Two waves of distortion
end

if in.FitPhi4
    OptionsExp      = [OptionsExp, ',''Phi4'',Phi4/',num2str(in.ScalePhi4)];      
    if isnan(in.Phi40)
        StartPoint  = [StartPoint,  0*in.ScalePhi4];
    else
        StartPoint  = [StartPoint,  in.Phi40*in.ScalePhi4];   
    end
    Lower           = [Lower,   -4*pi*in.ScalePhi4];
    Upper           = [Upper,    4*pi*in.ScalePhi4];                       % Two waves of distortion
end

if in.FitPhi5
    OptionsExp      = [OptionsExp, ',''Phi5'',Phi5/',num2str(in.ScalePhi5)];      
    if isnan(in.Phi50)
        StartPoint  = [StartPoint,  0*in.ScalePhi5];
    else
        StartPoint  = [StartPoint,  in.Phi50*in.ScalePhi5];   
    end
    Lower           = [Lower,   -4*pi*in.ScalePhi5];
    Upper           = [Upper,    4*pi*in.ScalePhi5];                       % Two waves of distortion
end

if in.FitPhi6
    OptionsExp      = [OptionsExp, ',''Phi6'',Phi6/',num2str(in.ScalePhi6)];      
    if isnan(in.Phi60)
        StartPoint  = [StartPoint,  0*in.ScalePhi6];
    else
        StartPoint  = [StartPoint,  in.Phi60*in.ScalePhi6];   
    end
    Lower           = [Lower,   -4*pi*in.ScalePhi6];
    Upper           = [Upper,    4*pi*in.ScalePhi6];                       % Two waves of distortion
end

if in.FitThetaOffset
    ThetaOffsetExp  = ['x - ThetaOffset/',num2str(in.ScaleThetaOffset)];   % Include x(theta) offset
    if isnan(in.ThetaOffset0)
        StartPoint  = [StartPoint,  0*in.ScaleThetaOffset];
    else
        StartPoint  = [StartPoint, in.ThetaOffset0*in.ScaleThetaOffset];
    end
    Lower           = [Lower,       2*in.ScaleThetaOffset];
    Upper           = [Upper,       2*in.ScaleThetaOffset];                 % +/- 2 degrees
end

if in.Fitalpha      % Done with upper case, so alpha is next
    OptionsExp      = [OptionsExp, ',''alpha'',alpha/',num2str(in.Scalealpha)];  
    if isnan(in.alpha0)
        StartPoint  = [StartPoint,  0*in.Scalealpha];
    else
        StartPoint  = [StartPoint,  in.alpha0*in.Scalealpha];   
    end
    Lower           = [Lower,       0*in.Scalealpha];
    Upper           = [Upper,    10^5*in.Scalealpha];                      % 1/e decay in 10 microns
end

if in.FitdPhi2
    OptionsExp      = [OptionsExp, ',''dPhi2'',dPhi2/',num2str(in.ScaledPhi2)];      
    if isnan(in.dPhi20)
        StartPoint  = [StartPoint,  0*in.ScaledPhi2];
    else
        StartPoint  = [StartPoint,  in.dPhi20*in.ScaledPhi2];   
    end
    Lower           = [Lower,   -4*pi*in.ScaledPhi2];
    Upper           = [Upper,    4*pi*in.ScaledPhi2];                      % Two waves of distortion
end

if in.FitdPhi3
    OptionsExp      = [OptionsExp, ',''dPhi3'',dPhi3/',num2str(in.ScaledPhi3)];
    if isnan(in.dPhi30)
        StartPoint  = [StartPoint,  0*in.ScaledPhi3];
    else
        StartPoint  = [StartPoint,  in.dPhi30*in.ScaledPhi3];
    end
    Lower           = [Lower,   -2*pi*in.ScaledPhi3];
    Upper           = [Upper,    2*pi*in.ScaledPhi3];                      % One wave of distortion
end

if in.FitdPhi4
    OptionsExp      = [OptionsExp, ',''dPhi4'',dPhi4/',num2str(in.ScaledPhi4)];
    if isnan(in.dPhi40)
        StartPoint  = [StartPoint,  0*in.ScaledPhi4];
    else
        StartPoint  = [StartPoint,  in.dPhi40*in.ScaledPhi4];
    end
    Lower           = [Lower,   -1*pi*in.ScaledPhi4];
    Upper           = [Upper,    1*pi*in.ScaledPhi4];                      % Half wave of distortion
end

if in.Fitn1
    OptionsExp      = [OptionsExp, ',''n1'',n1/',num2str(in.Scalen1)];
    if isnan(in.n10)
        StartPoint  = [StartPoint,  0.001*in.Scalen1];                     % Start at n1 = 10^-3
    else
        StartPoint  = [StartPoint,  in.n10*in.Scalen1];
    end
    Lower           = [Lower,       0*in.Scalen1];
    Upper           = [Upper,     0.1*in.Scalen1];                         % Max of n1 = 0.1
end

if in.FitdPhi4
    OptionsExp      = [OptionsExp, ',''z04OverL'',z04OverL/',num2str(in.Scalez04OverL)];
    if isnan(in.z04OverL0)
        StartPoint      = [StartPoint,.33*in.Scalez04OverL];             % Nulls are equally spaced 
    else
        StartPoint  = [StartPoint,in.z04OverL0*in.Scalez04OverL];       
    end
    Lower           = [Lower,       0*in.Scalez04OverL];
    Upper           = [Upper,     0.5*in.Scalez04OverL];                 % Max 1/2, so between 0 and L/2
end

%--------------------------------------------------------------------------
% Create pass thru parameter strings.  These are all the controls for the
% BraggTransmission function that *aren't* being fit.
%--------------------------------------------------------------------------
ParamExp = '';
for ip = 1:length(in.PassThruParams)
    ParamExp = [ParamExp,'''', in.PassThruParams{ip}, ''',', num2str(in.(in.PassThruParams{ip}))];
    if ip ~= length(in.PassThruParams)
        ParamExp = [ParamExp,','];        % Add comma for all but last entry
    end
end

% Same but for fit coefficients that have been turned off but starting value supplied.
CoefExp = '';
for ip = 1:length(in.PassThruCoeffs)
    CoefExp = [CoefExp,'''', in.PassThruCoeffs{ip}, ''',', num2str(in.([in.PassThruCoeffs{ip},'0']))];
    if ip ~= length(in.PassThruCoeffs)
        CoefExp = [CoefExp,','];        % Add comma for all but last entry
    end
end

%--------------------------------------------------------------------------
% Create fit type expression string and call the fittype matlab function
%--------------------------------------------------------------------------
FitTypeExp= [EffOffsetExp, 'BraggTransmission(', ThetaOffsetExp,',',ParamExp]; % Create initial ft expression
if ~isempty(CoefExp), FitTypeExp= [FitTypeExp,',',CoefExp]; end        % Add coefficient pass through

if ~isempty(OptionsExp)
    FitTypeExp = [FitTypeExp,OptionsExp];               % Add options controlling fit parameters
end
FitTypeExp = [FitTypeExp,')'];                          % Add closing paren

out.ft	= fittype(FitTypeExp, 'independent', 'x', 'dependent', 'y' );

%--------------------------------------------------------------------------
% Initialize fit options structure
%--------------------------------------------------------------------------
opts                = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display        = 'off';
opts.StartPoint     = StartPoint;
opts.Lower          = Lower;
opts.Upper          = Upper;

%--------------------------------------------------------------------------
% Fit model to data.
%--------------------------------------------------------------------------
[xData, yData]              = prepareCurveData( DeltaThetaExt, EffVsDeltaTheta );    % Clean up fit inputs
[out.fitresult, out.gof]    = fit( xData, yData, out.ft, opts );

%--------------------------------------------------------------------------
% Extract results and confidence intervals and store in "out"
%--------------------------------------------------------------------------
CoefficientNames    = coeffnames(out.ft);           % Get names of fit coefficients in order
CoefficientConfInt  = confint(out.fitresult,0.95);  % Confidence intervals in same order

for ic = 1:length(CoefficientNames)
    Value = out.fitresult.(CoefficientNames{ic});                       % Extract value
    out.(CoefficientNames{ic}) = Value/in.(['Scale',CoefficientNames{ic}]); % Create field in out with this name and value
    out.([CoefficientNames{ic},'UpperCI']) = CoefficientConfInt(2,ic)/in.(['Scale',CoefficientNames{ic}]);  % Create confidence interval field
end

end % FitExpToTheory

%==========================================================================
% Plot results
%==========================================================================
function PlotResults(DeltaThetaExt,EffVsDeltaTheta,in,out)

figure;hold off;
plot(DeltaThetaExt,EffVsDeltaTheta,'k.');        % Data
hold on;
plot(out.fitresult,'r');                            % Fit
legend('hide');         % Supress "fitted curve" legend

% set(gcf,'color','w');
% set(gca,'FontSize',12);
% set(gca,'FontName','Arial');
xlabel('\Delta{\theta}_{0} [^o]');
ylabel('\eta');
title(in.PlotTitle);

end % PlotResults

%==========================================================================
% Plot fit results on current figure.  Build string with newlines, then
% place on figure.
%==========================================================================
function PlotFitResults(in,out)

% Where to locate table
xlimits     = xlim;
ylimits     = ylim;
EffLoc      = 0.95*ylimits(2);   % 1.05        
ThetaLoc    = 0.12*xlimits(2);    % 0.5

% Fit quality
ResultsStr = ['R^2        = ',num2str(round(out.gof.rsquare,3))];

% Born or not? % Commented out by AS 2025-03-03
% ResultsStr = [ResultsStr,newline,'Born       = ',num2str(in.Born)];

% Thickness
ResultsStr = [ResultsStr,newline,...
    'L / {\mu}m   = ',num2str(round(out.L/1E-6,2)),...
    ' {\pm} ',num2str(round((out.LUpperCI-out.L)/1E-6,2))];

% Delta n
ResultsStr = [ResultsStr,newline,...
    'n_1 / 10^{-3} = ',num2str(round(out.n1/0.001,2)),...
    ' {\pm} ',num2str(round((out.n1UpperCI-out.n1)/0.001,2))];

% Average delta n
if in.FitA1 || in.FitA2  || in.FitA3  || in.FitA4  || in.FitA5 || in.FitA6 
    if isfield(out, 'A1'),     A1    = out.A1;     else;  A1    = 0; end        % Values not defined if not fit
    if isfield(out, 'A3'),     A3    = out.A3;     else;  A3    = 0; end        % Even terms have average of zero
    if isfield(out, 'A5'),     A5    = out.A5;     else;  A5    = 0; end
    Abar = 1 + A1*2/pi + A3*2/3/pi + A5*2/5/pi;                                 % Average amplitude in z

    ResultsStr = [ResultsStr,newline,...    
        '<n_1>/10^{-3} = ',num2str(round(out.n1*Abar/0.001,2)),...
        ' {\pm} ',num2str(round((out.n1UpperCI-out.n1)*Abar/0.001,2))];
end

% Optional fit coefficients
if in.FitdPhi2
    ResultsStr = [ResultsStr,newline,...
        '{\delta}{\phi}_2 / 2{\pi}  = ',num2str(round(out.dPhi2/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.dPhi2UpperCI-out.dPhi2)/(2*pi),2))];
end

if in.FitdPhi3
    ResultsStr = [ResultsStr,newline,...
        '{\delta}{\phi}_3 / 2{\pi}  = ',num2str(round(out.dPhi3/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.dPhi3UpperCI-out.dPhi3)/(2*pi),2))];
end

if in.FitdPhi4
    ResultsStr = [ResultsStr,newline,...
        '{\delta}{\phi}_4 / 2{\pi}  = ',num2str(round(out.dPhi4/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.dPhi4UpperCI-out.dPhi4)/(2*pi),2))];
    ResultsStr = [ResultsStr,newline,...
        'z_{04} / L    = ',num2str(round(out.z04OverL,2)),...
        ' {\pm} ',num2str(round((out.z04OverLUpperCI-out.z04OverL),2))];
end

if in.FitPhi1
    ResultsStr = [ResultsStr,newline,...
        '{\phi}_1 / 2{\pi}   = ',num2str(round(out.Phi1/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.Phi1UpperCI-out.Phi1)/(2*pi),2))];
end

if in.FitPhi2
    ResultsStr = [ResultsStr,newline,...
        '{\phi}_2 / 2{\pi}   = ',num2str(round(out.Phi2/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.Phi2UpperCI-out.Phi2)/(2*pi),2))];
end

if in.FitPhi3
    ResultsStr = [ResultsStr,newline,...
        '{\phi}_3 / 2{\pi}   = ',num2str(round(out.Phi3/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.Phi3UpperCI-out.Phi3)/(2*pi),2))];
end

if in.FitPhi4
    ResultsStr = [ResultsStr,newline,...
        '{\phi}_4 / 2{\pi}   = ',num2str(round(out.Phi4/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.Phi4UpperCI-out.Phi4)/(2*pi),2))];
end

if in.FitPhi5
    ResultsStr = [ResultsStr,newline,...
        '{\phi}_5 / 2{\pi}   = ',num2str(round(out.Phi5/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.Phi5UpperCI-out.Phi5)/(2*pi),2))];
end

if in.FitPhi6
    ResultsStr = [ResultsStr,newline,...
        '{\phi}_6 / 2{\pi}   = ',num2str(round(out.Phi6/(2*pi),2)),...
        ' {\pm} ',num2str(round((out.Phi6UpperCI-out.Phi6)/(2*pi),2))];
end

if in.Fitalpha
    ResultsStr = [ResultsStr,newline,...
        '{\alpha} L         = ',num2str(round(out.alpha*out.L,2)),...
        ' {\pm} ',num2str(round((out.alphaUpperCI-out.alpha)*out.L,2))];        
end

if in.FitA1
    ResultsStr = [ResultsStr,newline,...
        'A_1          = ',num2str(round(out.A1,2)),...
        ' {\pm} ',num2str(round(out.A1UpperCI-out.A1,2))];
end

if in.FitA2
    ResultsStr = [ResultsStr,newline,...
        'A_2          = ',num2str(round(out.A2,2)),...
        ' {\pm} ',num2str(round(out.A2UpperCI-out.A2,2))];
end

if in.FitA3
    ResultsStr = [ResultsStr,newline,...
        'A_3          = ',num2str(round(out.A3,2)),...
        ' {\pm} ',num2str(round(out.A3UpperCI-out.A3,2))];
end

if in.FitA4
    ResultsStr = [ResultsStr,newline,...
        'A_4          = ',num2str(round(out.A4,2)),...
        ' {\pm} ',num2str(round(out.A4UpperCI-out.A4,2))];
end

if in.FitA5
    ResultsStr = [ResultsStr,newline,...
        'A_5          = ',num2str(round(out.A5,2)),...
        ' {\pm} ',num2str(round(out.A5UpperCI-out.A5,2))];
end

if in.FitA6
    ResultsStr = [ResultsStr,newline,...
        'A_6          = ',num2str(round(out.A6,2)),...
        ' {\pm} ',num2str(round(out.A6UpperCI-out.A6,2))];
end

if in.FitEffOffset
    ResultsStr = [ResultsStr,newline,...
        '{\Delta}{\eta}          = ',num2str(round(out.EffOffset,2)),...
        ' {\pm} ',num2str(round((out.EffOffsetUpperCI-out.EffOffset),2))];
end

if in.FitThetaOffset
    ResultsStr = [ResultsStr,newline,...
        '{\Delta}{\theta}       = ',num2str(round(out.ThetaOffset,2)),...
        ' {\pm} ',num2str(round((out.ThetaOffsetUpperCI-out.ThetaOffset),2))];
end

% Show on plot
text(ThetaLoc,EffLoc,ResultsStr,'FontSize',10,'FontName','Arial',...
    'EdgeColor','k','BackgroundColor','w','VerticalAlignment','top');

end % PlotFitResults

%==========================================================================
% Plot disortion on current figure. 
%==========================================================================
function PlotPhi(~,out)

% Get fit parameters
if isfield(out,'dPhi2'),    dPhi2    = out.dPhi2;    else; dPhi2    = 0; end
if isfield(out,'dPhi3'),    dPhi3    = out.dPhi3;    else; dPhi3    = 0; end
if isfield(out,'dPhi4'),    dPhi4    = out.dPhi4;    else; dPhi4    = 0; end
if isfield(out,'z04OverL'), z04OverL = out.z04OverL; else; z04OverL = 0; end

if isfield(out, 'Phi1'),     Phi1    = out.Phi1;     else;  Phi1    = 0; end
if isfield(out, 'Phi2'),     Phi2    = out.Phi2;     else;  Phi2    = 0; end
if isfield(out, 'Phi3'),     Phi3    = out.Phi3;     else;  Phi3    = 0; end
if isfield(out, 'Phi4'),     Phi4    = out.Phi4;     else;  Phi4    = 0; end
if isfield(out, 'Phi5'),     Phi5    = out.Phi5;     else;  Phi5    = 0; end
if isfield(out, 'Phi6'),     Phi6    = out.Phi6;     else;  Phi6    = 0; end

% Calculate phase profile vs depth
zOverL  = 0:.02:1;
if Phi1  == 0 && Phi2  == 0 && Phi3  == 0 && Phi4  == 0 && Phi5  == 0 && Phi6  == 0
    Phi = PhiPoly(zOverL,dPhi2,dPhi3,dPhi4,z04OverL);
else
    Phi = FourierSeries(zOverL,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6);
end

% Create new axes
axes('Position',[.22 .7 .2 .2])
box on

plot(zOverL,Phi/(2*pi));
set(gcf,'color','w');
set(gca,'FontSize',12);
set(gca,'FontName','Arial');
xlabel('z / L');
%title('Distortion');
ylabel('\phi / 2{\pi}');

end % PlotPhi

%==========================================================================
% Plot attenuation on current figure. 
%==========================================================================
function PlotA(~,out)

% Get fit parameters
if isfield(out, 'A1'),     A1    = out.A1;     else;  A1    = 0; end
if isfield(out, 'A2'),     A2    = out.A2;     else;  A2    = 0; end
if isfield(out, 'A3'),     A3    = out.A3;     else;  A3    = 0; end
if isfield(out, 'A4'),     A4    = out.A4;     else;  A4    = 0; end
if isfield(out, 'A5'),     A5    = out.A5;     else;  A5    = 0; end
if isfield(out, 'A6'),     A6    = out.A6;     else;  A6    = 0; end

% Calculate phase profile vs depth
zOverL  = 0:.02:1;

% Note addition of fixed "1" term.  See BraggTransmission/DETransmissionBornFourierDistortedAndAttenuatedNumeric
A = 1 + FourierSeries(zOverL,A1,A2,A3,A4,A5,A6);  

% Create new axes
axes('Position',[.22 .45 .2 .2])
box on

plot(zOverL,A);
set(gcf,'color','w');
set(gca,'FontSize',12);
set(gca,'FontName','Arial');
xlabel('z / L');
%title('Attenuation');
ylabel('A');

end % PlotA

%==========================================================================
% Polynomial distortion of phase.  Copied from BraggTransmission since
% awkward to pass back through fit.
%==========================================================================
function Phi = PhiPoly(zOverL,dPhi2,dPhi3,dPhi4,z04OverL)

% Normalization of quartic depends on z04OverL since two possible extrema
% for a 4th order polynomial
Norm4 = max([(1-2*z04OverL)^2/16,(z04OverL-1)^2 * z04OverL^2 / 4 ]);

Phi = dPhi2 * 4 *          zOverL.*(1-zOverL) + ...
      dPhi3 * 12*sqrt(3) * zOverL.*(zOverL-.5).*(zOverL-1) + ...
      dPhi4 *              zOverL.*(zOverL-z04OverL).*(zOverL+z04OverL-1).*(zOverL-1)/Norm4;

end % PhiPoly

%==========================================================================
% Fourier series distortion of phase or amplitude
%==========================================================================
function Phi = FourierSeries(zOverL,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6)

Phi = Phi1 * sin(2*pi*zOverL * 1 /2) + ...
      Phi2 * sin(2*pi*zOverL * 2 /2) + ...
      Phi3 * sin(2*pi*zOverL * 3 /2) + ...
      Phi4 * sin(2*pi*zOverL * 4 /2) + ...
      Phi5 * sin(2*pi*zOverL * 5 /2) + ...
      Phi6 * sin(2*pi*zOverL * 6 /2);

end % PhiFourier

%==========================================================================
% Plot k space of writing and reading.  All quantities rendered internal to
% the material.
%==========================================================================
function PlotKSpace(in,out)

figure;

HolographicRenderer2D('RefWrtTheta0',in.RefWrtTheta*pi/180,'ObjWrtTheta0',in.ObjWrtTheta*pi/180,...
    'lambdaWrt0',in.lambdaWrt0,'lambdaRed0',in.lambdaRed0,'nWrt0',in.nWrt,'nRed0',in.nRed,...
    'L',out.L,'H',5*out.L);

HolographicRenderer2D('RefWrtTheta0',in.RefWrtTheta*pi/180,'ObjWrtTheta0',in.ObjWrtTheta*pi/180,...
    'lambdaWrt0',in.lambdaWrt0,'lambdaRed0',in.lambdaRed0,'nWrt0',in.nWrt,'nRed0',in.nRed,...
    'L',out.L,'H',5*out.L,...
    'Reading',true,'DrawHorizAxis',false,'DrawVertAxis',false,'DrawUncertainty',true);

end % PlotKSpace
