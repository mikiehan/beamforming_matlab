%% Examples on how to calculate array factor / beampattern (Delay and sum algorithm)
close all
clc
clear

nrow = 1;
ncol = 3;
shiftn = 0;
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.5]);
title('normalized beampattern')

%% parameters
% array type
arraytype = 'uca'; %'ula';  % 'ula', 'uca', 'nla'

% Wave-frequency and wave-speed
f = [1e3 3e3 8e3 19e3];
c = 340;

% Scanning angles
thetaScanningAngles = 0:360;
phiScanningAngles = 90;

% Steering angle
thetaSteeringAngle = 90;
phiSteeringAngle = 0;

%% Create mic array
% Position of sensors and weighting of mic array
if strcmp(arraytype, 'nla')
    % Create a linear array
    radius = 0.01;
    xPos = radius * cumsum([0, 3, 2, 3, 2, 3, 2, 3]);
    yPos = zeros(1, numel(xPos));
    zPos = zeros(1, numel(xPos));
    nElements = numel(xPos);   
elseif strcmp(arraytype, 'ula')
    % Create a linear array
    radius = 0.01;
    xPos = radius * cumsum([0, 14/10*ones(1,5)]);
    yPos = zeros(1, numel(xPos));
    zPos = zeros(1, numel(xPos));
    nElements = numel(xPos);   
elseif strcmp(arraytype, 'uca')
    % Create circular array
    nElements = 8;
    radius = 0.0383;
    [xPos, yPos] = pol2cart(linspace(0,2*pi-2*pi/nElements,nElements),ones(1,nElements)*radius);
    zPos = zeros(1, numel(xPos));
end
elementWeights = ones(1, numel(xPos))/sqrt(numel(xPos));

%% Plot array geometry
%Calculate and plot the array pattern for various frequencies
axGeometry = subplot(nrow, ncol, 1+shiftn);
scatter(axGeometry, xPos, yPos, 'filled')
title(axGeometry, 'Array geometry','FontWeight','Normal')
% axis(axGeometry, [-radius-0.1 radius+0.1 -radius-0.1 radius+0.1])
axis(axGeometry, 'square')
grid(axGeometry, 'minor')
xlabel('m')
ylabel('m')

%% Plot the beampattern for various frequencies with plotBeampattern()
subplot(nrow, ncol, [1,3])
plotBeampattern(xPos, yPos, zPos, elementWeights, f, c, thetaSteeringAngle, phiScanningAngles, 50, 'rect',nrow, ncol,shiftn)

%% Plot array factor from 0 to 20 kHz
% Preallocating for speed
ff = 0:10:20e3;
W_all = zeros(numel(ff), numel(thetaScanningAngles));
% Calculate array factor
for k = 1:length(ff)
    ei = steeringVector(xPos, yPos, zPos, ff(k), c, thetaScanningAngles, phiScanningAngles);
    AF = arrayFactor(elementWeights, ei);
    W_all(k,:) = 20*log10(AF);
end

% Don't display values below a certain threshold
dynRange = 30;
W_all(W_all<-dynRange) = NaN;

% Plot array factor
ax = subplot(nrow, ncol, 3+shiftn);
imagesc(ax, ff*1e-3,thetaScanningAngles, W_all')
xlabel(ax, 'kHz')
ylabel(ax, '\theta')
colorbar(ax, 'northoutside', 'direction', 'reverse')

