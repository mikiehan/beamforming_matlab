amp = 1; % if greater than 1 then data will be clipped so keep 1 or under
duration = 1; % 10 seconds
fc = 20e3; % Hz
sampling_freq = fc*2; % sampling frequency
t = 0:1/sampling_freq:duration;
t = t(1:end-1)';
txSignal = amp*sin(2*pi*fc*t);
s_Pos = [0, 5, 0]; % source position in meter
c = 340;
%c = physconst('LightSpeed');
radius = 0.05; 
mic = phased.OmnidirectionalMicrophoneElement(...
     'FrequencyRange',[20 fc+5e3]);
P = 10;
incidentAngle = [90;0];
wavelength = c/fc;

array = phased.UCA('NumElements', P ,'Radius', radius ,'Element',mic);
pos = getElementPosition(array);

m_xPos = pos(1, :); % x-corrd
m_yPos = pos(2, :); % y-coord
m_zPos = pos(3, :); % z-coord 

% m_xPos = [0, -38.13, -20.98, 11.97, 35.91, 32.81, 5.00, -26.57] ./ 1000 * 2;
% m_yPos = [0, 3.58, 32.04, 36.38, 13.32, -19.77, -37.97, -27.58] ./ 1000 * 2;
% m_zPos = zeros(1, P);
% P = length(m_xPos);
%[m_gamma, m_l] = cart2pol(m_xPos, m_yPos);


% distance = zeros(1, P);
% r_amp = zeros(1, P);
% for i=1:P
%     X = [s_Pos(1), s_Pos(2); m_xPos(i), m_yPos(i)];
%     distance(i) = pdist(X, 'euclidean');
%     r_amp(i) = (wavelength/(4*pi*distance(i)))^2;
% end
% 
% rxSignal = txSignal * r_amp;
rxSignal = collectPlaneWave(array,txSignal,incidentAngle,fc,c);


noise = 0.0000000001*(randn(size(rxSignal)) + 1j*randn(size(rxSignal)));
rxSignal = rxSignal + noise; % adding some noise

beamformer = phased.PhaseShiftBeamformer('SensorArray',array,...
    'PropagationSpeed',c,'OperatingFrequency',fc,...
    'Direction',incidentAngle,'WeightsOutputPort',true);

% steervec = phased.SteeringVector('SensorArray',array,...
%     'PropagationSpeed',c);
% beamformer = phased.LCMVBeamformer('Constraint',steervec(fc,incidentAngle),'DesiredResponse',1,'WeightsOutputPort',true);

[y,w] = beamformer(rxSignal);

figure;
plot(t,real(rxSignal(:,1)),'r:',t,real(y),'b')
xlim([0.01 0.03]);
xlabel('Time')
ylabel('Amplitude')
legend('Original','Beamformed')
% 
% ebi_full = phased.SteeringVector('SensorArray',array,...
%     'PropagationSpeed',c);
% %ebi = ebi_full(fc, phitheta2azel([0;90]));
% ebi = ebi_full(fc, phitheta2azel([0;180]));
% % ebi = ebi_full(fc, phitheta2azel([0;0]));
% ebi = reshape(ebi, 1, 1, 8);
% 
% 
% rxSignal = rxSignal';
% wi = weightingVectorMVDR(rxSignal, ebi);
% beamformed_signal = zeros(1, 1, length(rxSignal));
% for theta_angle = 1:1
%     for phi_angle = 1:1
%         weights = zeros(1, P);
%         for channel = 1:P
%             weights(channel) = wi(theta_angle,phi_angle,channel);
%         end
%         beamformed_signal(theta_angle, phi_angle,:) = weights *  rxSignal;
%     end
% end
% 
% beamformed_at_aoa = real(beamformed_signal);
% beamformed_at_aoa = reshape(beamformed_at_aoa, [], 1);
% t=t';
% figure;
% plot(t,real(rxSignal(3,:)),'r:',t,real(beamformed_at_aoa),'b')
% xlim([0 0.001])
% xlabel('Time')
% ylabel('Amplitude')
% legend('Original','Beamformed')

