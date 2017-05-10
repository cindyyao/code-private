function PhotonCatchRate = PhotonCatchRateLED(Powers, radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START HERE

% Set physical constants
SpeedOfLight = 3e8;
PlanksConstant = 6.626e-34;
WaveLengthRange = 778:-2:380; % in nm
OpsinPhotosensitivity = 9.6e-21; % m^2 
wavelengths = 400:1:800;

%% UDT
%UDT Information
cd /Volumes/lab/Experiments/Calibration/2013-07-19/ 
load cal_udt_spectrum;



%% LED
% Powers = 4.5e-12; % Watts: Power made with 
spot_area = radius^2 * pi; % radius = 0.0025
cd /Volumes/lab/Experiments/Calibration/LED_Spectra/ 
load LED_spectrum
LED_spectrum = LED_spectrum ./ norm(LED_spectrum);

RGBMatrix = LED_spectrum;


%% LED calibration

% Normalize and scale the RGBMAtrix
TruePowerScalers = Powers ./ (cal_udt_spectrum(:,2)' * RGBMatrix');
CalibratedRGBMatrix = RGBMatrix' * TruePowerScalers;

SummedMonitorSpectra = CalibratedRGBMatrix'; % this is to make it compatible with a multi-primary source

%Sanity Check: CheckPower should equal the sum of Powers
Powers
CheckPower = dot(CalibratedRGBMatrix, repmat(cal_udt_spectrum(:,2), 1, length(Powers)))


%% Rod calculation begins here:
% Rod spectral sensitivty fit with Baylor nomogram
WL = 800:-1:400;
WL = WL.*0.001;
a0 = -5.2734;
a1 = -87.403;
a2 = 1228.4;
a3 = -3346.3;
a4 = -5070.3;
a5 = 30881;
a6 = -31607;
lambdaMax = 491;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
RodPhotonSensitivity = 10.^LogPhotonSensitivity;

WavelengthsCorrected = wavelengths * 1e-9; % nm -> m
Intensity = (SummedMonitorSpectra .* repmat(WavelengthsCorrected, length(Powers), 1)) ./ (PlanksConstant * SpeedOfLight);
PhotonFlux = Intensity ./ spot_area;
AbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux',repmat(RodPhotonSensitivity([401:-1:1])', 1, length(Powers)));

%RodCollectingArea = 1.2e-12; % MONKEY (Baylor 1984, Table 1) m^2
RodCollectingArea = 0.5e-12; % MOUSE m^2
EffectivePhotonFlux = dot(PhotonFlux',repmat(RodPhotonSensitivity([401:-1:1])', 1, length(Powers)));
PhotonCatchRate = EffectivePhotonFlux * RodCollectingArea;



