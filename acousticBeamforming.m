Nr = 8;        %# of antenna element of Rx
fc = 20e3;   %central carrier frequency
c = 340;
lambda = c/fc;  %wavelength


radius = 0.05; % 3.83 cm
rxarray = phased.UCA('NumElements', Nr ,'Radius', radius);

s_Pos = [0, 1, 0]; % source position in meter

pos = getElementPosition(rxarray);
m_xPos = pos(1, :); % x-coord
m_yPos = pos(2, :); % y-coord
m_zPos = pos(3, :); % z-coord

Elem = phased.OmnidirectionalMicrophoneElement;
Elem.FrequencyRange = [0 fc];
rxarray.Element = Elem;

incidentAngle = [90; 0]
ang_azi_rx = [90];
SteeringAngles_Rx = [ang_azi_rx;zeros(1,length(ang_azi_rx))];
PhaseShiftBits_Rx = 0;

duration = 1;
sampling_freq = fc*2; % sampling frequency
t = 0:1/sampling_freq:duration;
t = t(1:end-1)';

distance = zeros(1, Nr);
r_attn = zeros(1, Nr); % attenuation
r_delay = zeros(1, Nr); % delay
rx = zeros(length(t), Nr);
for i=1:Nr
    X = [s_Pos(1), s_Pos(2); m_xPos(i), m_yPos(i)];
    distance(i) = pdist(X, 'euclidean');
    r_delay(i) = distance(i)/c;
    r_attn(i) = (1/(distance(i)))^2;
    rx(:, i) = r_attn(i)*sin(2*pi*fc*t + r_delay(i));
end

noise = 0; % 0.1*(randn(size(rx)) + 1j*randn(size(rx)));
rxSignal = rx + noise; % adding some noise


% Calculate Rx Steering Weights
wr = zeros(getNumElements(rxarray), size(SteeringAngles_Rx, 2));
for i = 1:size(SteeringAngles_Rx, 2)
    steervec = phased.SteeringVector('SensorArray', rxarray,...
        'PropagationSpeed', c, ...
        'NumPhaseShifterBits', PhaseShiftBits_Rx);
    beamformer = phased.LCMVBeamformer('Constraint',steervec(fc,incidentAngle),'DesiredResponse',1,'WeightsOutputPort',true);
    [y,w] = beamformer(rxSignal);
    
    ebi = steervec(fc, incidentAngle);
    ebi = reshape(ebi, 1, 1, 8);
    
    wi = weightingVectorMVDR(rxSignal', ebi);
    y2 = rxSignal*squeeze(wi);
    
%     beamformed_at_aoa = real(beamformed_signal);
%     beamformed_at_aoa = reshape(beamformed_at_aoa, [], 1);
%     
    figure;
    plot(t,real(rxSignal(:,2)),'b:',t,real(y),'r', t,real(rx(:,2)),'c', t, real(y2), 'b');
    xlim([0 0.001]);
    xlabel('Time')
    ylabel('Amplitude')
    legend('Raw','Beamformed', 'Original no noise', 'Beamformed 2')
end

figure;
title('Rx Beam Pattern (our)')
format = 'polar';

cutAngle = 0;
pattern(rxarray, fc, [-180:-90, 90:180], cutAngle, 'PropagationSpeed', c,...
    'Type', 'efield', 'CoordinateSystem', format ,'weights',squeeze(wi));

figure;
title('Rx Beam Pattern (matlab)')
format = 'polar';

cutAngle = 0;
pattern(rxarray, fc, [-180:-90, 90:180], cutAngle, 'PropagationSpeed', c,...
    'Type', 'efield', 'CoordinateSystem', format ,'weights',w);

