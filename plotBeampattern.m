function [] = plotBeampattern(xPos, yPos, zPos, weights, f, c, thetaSteerAngle, phiScanAngle, dynRange, plotType, nrow,ncol,shiftn)
%plotBeampattern - plots the beampattern for various frequencies
%
%plotBeampattern(xPos, yPos, zPos, weights, f, c, thetaSteerAngle, phiScanAngle, dynRange, plotType)
%
%IN
%xPos             - 1xP vector of x-positions [m]
%yPos             - 1xP vector of y-positions [m]
%zPos             - 1xP vector of z-positions [m]
%weights          - 1xP vector of element weights (optional, default uniform weighting)
%f                - Wave frequency [Hz] (optional, default 0.5, 1, 1.5, 3 kHz)
%c                - Speed of sound [m/s] (optional, default 340 m/s)
%thetaSteerAngle  - 1x1 theta steering angle [degrees]  (optional)
%phiScanAngle     - Angle slice to show, 0 for xz and 90 for yz view  (optional)
%dynRange         - Dynamic range in plot [dB]  (optional)
%plotType         - Use 'rect' or 'polar' (optional)
%
%OUT
%[]               - The figure plot
%

if ~exist('plotType','var')
    plotType = 'full';
end

if ~exist('dynRange','var')
    dynRange = 50;
end

if ~exist('phiScanAngle', 'var')
    phiScanAngle = 90;
end

if ~exist('thetaSteerAngle', 'var')
    thetaSteerAngle = 0;
end

if ~exist('c', 'var')
    c = 340;
end

if ~exist('f', 'var')
    f = [0.5 1 1.5 3]*1e3;
end

if ~exist('weights', 'var')
    weights = ones(1, numel(xPos))/numel(xPos);
end

if ~exist('shiftn', 'var')
    shiftn = 0;
end

maxweight = 20*log10(sum(weights));
scaling = 1+maxweight/dynRange;
dBTicks = [maxweight:-20:0 0 -(20:20:dynRange-10)];
angleTicks = 0:30:180;

%Scanning angles
thetaScanAngles = 0:0.01:360;

%Linewidth in plot
lwidth = 1.5;

% Plot beampattern
% bpFig = figure;
% bpFig.Color = [1 1 1];

%Rectangular plot axis
if plotType ~= "polar"
    axRectangular = subplot(nrow, ncol, 1+shiftn);
    hold(axRectangular, 'on')
    xlabel(axRectangular, 'Angle (deg)')
    ylabel(axRectangular, 'Attenuation (dB)')
    grid(axRectangular, 'on')
    axis(axRectangular, [thetaScanAngles(1) thetaScanAngles(end) -dynRange maxweight])
    % axRectangular.YTick = [-50 -45 -40 -35 -30 -25 -20 -15 -10 -6 -3 0];
    % axRectangular.XTick = [-90 -80 -70 -60 -50 -40 -30 -20 -10 0 10 20 30 40 50 60 70 80 90];
end

if plotType ~= "rect"
    %Polar plot axis
    axPolar = subplot(nrow, ncol, 2+shiftn);
    axPolar.Visible = 'off';
    hold(axPolar, 'on')
    axis(axPolar, 'equal')
    ylim(axPolar, [0 dynRange]*scaling)
    xlim(axPolar, [-dynRange dynRange]*scaling)
    
    %Set background color for polar plot to white (inside half circle)
    patch('XData', cos(0:pi/50:2*pi) * dynRange * scaling, ...
        'YData', sin(0:pi/50:2*pi) * dynRange * scaling,...
        'FaceColor', [1 1 1], ...
        'Parent', axPolar);

    %Plot angle spokes in polar plot
    for tick = angleTicks
        line(dynRange * scaling * [-sin(tick*pi/180) sin(tick*pi/180)], ...
            dynRange * scaling * [-cos(tick*pi/180) cos(tick*pi/180)], ...
            'LineStyle', '-', ...
            'Color', [1 1 1]*0.8, ...
            'LineWidth', 0.5, ...
            'Parent', axPolar);

        text((dynRange*scaling*1.08) * cos(tick*pi/180), ...
            (dynRange*scaling*1.08) * sin(tick*pi/180), ...
            [int2str(tick) '^\circ'],...
            'HorizontalAlignment', 'center', ...
            'fontSize', 10, ...
            'Parent', axPolar);
    end

    %Plot dB ticks in polar plot
    txtAngle = 10;
    for tick = dBTicks
        line(cos(0:pi/50:2*pi)*(dynRange+tick), sin(0:pi/50:2*pi)*(dynRange+tick), ...
            'LineStyle', '-', ...
            'Color', [1 1 1]*0.8, ...
            'LineWidth', 0.5, ...
            'Parent', axPolar);
    
    text((dynRange+tick)*cos(txtAngle*pi/180), ...
        (dynRange+tick)*sin(txtAngle*pi/180), ...
        ['  ' num2str(round(tick)) ' dB'], ...
        'fontsize',8, ...
        'Parent', axPolar);
    end
    line(cos(0:pi/50:2*pi)*dynRange*scaling, sin(0:pi/50:2*pi)*dynRange*scaling, ...
        'LineStyle', '-', ...
        'Color', [0 0 0], ...
        'LineWidth', 1, ...
        'Parent', axPolar);
    line([-dynRange*scaling dynRange*scaling], [0 0], ...
        'LineStyle', '-', ...
        'Color', [0 0 0], ...
        'LineWidth', 1, ...
        'Parent', axPolar);
end


%Calculate and plot the beampattern(s) in the figure
polarPlotHandles = [];
for ff = f
    ei = steeringVector(xPos, yPos, zPos, ff, c, thetaScanAngles, phiScanAngle);
    W = arrayFactor(weights, ei);
    W = 20*log10(W);
    W = reshape(W, 1, numel(W));
    
    % Rectangular plot
    if ff < 1e3
        displayName = [num2str(ff) ' Hz'];
    else
        displayName = [num2str(ff*1e-3) ' kHz'];
    end
    
    if plotType ~= "polar"
        plot(axRectangular, thetaScanAngles, W, 'linewidth', lwidth, 'DisplayName', displayName);
    end
    if plotType ~= "rect"
        % Polar plot
        xx = (W+dynRange) .* cos(thetaScanAngles*pi/180);
        yy = (W+dynRange) .* sin(thetaScanAngles*pi/180);
        p = plot(axPolar, xx, yy, 'linewidth', lwidth, 'DisplayName', displayName);
        polarPlotHandles = [polarPlotHandles p];
    end
end
if plotType ~= "polar"
	legend(axRectangular, 'show', 'Location','Eastoutside')
end
if plotType ~= "rect"
    legend(axPolar, polarPlotHandles, 'Location','Eastoutside')
end

% bpFig.Position = [500 200 540 600];

%Only show rectangular or polar plot if plotType is given
% switch plotType
%     case 'rect'
%         delete(axPolar)
%         % bpFig.Position = [500 200 540 300];
%         axRectangular.Position = [0.1300 0.1100 0.7750 0.7750];
%     case 'polar'
%         delete(axRectangular)
%         % bpFig.Position = [500 200 540 300];
%         axPolar.Position = [0.1300 0.1100 0.7750 0.7750];
% end
