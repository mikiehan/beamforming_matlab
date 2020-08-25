addpath('/Users/profhan/Research/fmcw_aoa_codes/matlab_codes/beamform_aoa_matlab/');
use_sim = false;
m_xPos = [0, -38.13, -20.98, 11.97, 35.91, 32.81, 5.00, -26.57] ./ 1000;
m_yPos = [0, 3.58, 32.04, 36.38, 13.32, -19.77, -37.97, -27.58] ./ 1000;
P = length(m_xPos);
m_zPos = zeros(1, P);
[m_gamma, m_l] = cart2pol(m_xPos, m_yPos);

s_xPos = 0;
s_yPos = 1;

fc = 19e3;
incidentAngle = [90;0];
c = 340; % (m/s) speed of sound

distance = zeros(1, P);
for i=1:P
    X = [s_xPos, s_yPos; m_xPos(i), m_yPos(i)];
    distance(i) = pdist(X, 'euclidean');
end
t = c/(max(distance) - min(distance));
disp(t);

mic = phased.OmnidirectionalMicrophoneElement(...
    'FrequencyRange',[20 20e3]);

nvec = [90, m_gamma(2:end)./pi * 180];
array = phased.ConformalArray(...
    'ElementPosition',[m_xPos;m_yPos;m_zPos],...
    'ElementNormal',[nvec;zeros(1,P)],...
    'Element', mic);
%array = phased.ULA('NumElements',5,'ElementSpacing',0.5);
%x = collectPlaneWave(array,xm,incidentAngle,fc,c);

figure;
viewArray(array, 'ShowNormals', true);
figure;
plotResponse(array, fc, c);

Fs = 48000; % sampling frequency
dt = 1/Fs;
duration = 0.5; % seconds
t=[0:dt:duration]';
t=t(1:end-1);

if (use_sim)
    amp = 1;
    signal=amp*sin(2*pi*fc*t);
    figure;
    plot(t(1:120), signal(1:120));
    
    figure;
    %%For one cycle get time period
    T = 1/fc ;
    % time step for one time period
    tt = 0:dt:T+dt ;
    d = sin(2*pi*fc*tt) ;
    plot(tt,d) ;
    
    % t = [0:.1:200]';
    % fr = .01;
    % signal = sin(2*pi*fr*t);
    
    x = collectPlaneWave(array,signal,incidentAngle,fc,c);
    noise = 0.001*(randn(size(x)) + 1j*randn(size(x)));
    rx = x + noise;
    
else
    inch = 48;
    trial = 1;
    prefix = sprintf('%s_%d%s_%d','fmcw_18k_20k_30', inch, 'inch_tgrhre', trial);
    %prefix = sprintf('%s_%d%s_%d', 'swadhin90',inch ,'inch_test', trial);
    display(prefix);
    [y0,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_0.wav'));
    [y1,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_1.wav'));
    [y2,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_2.wav'));
    [y3,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_3.wav'));
    [y4,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_4.wav'));
    [y5,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_5.wav'));
    [y6,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_6.wav'));
    [y7,Fs] = audioread(strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/' , prefix , '_channel_7.wav'));
    
    % apply FFT filter for each signal
    fmin = 18e3;
    fmax = 20e3;
    y0 = fftFilterSingleSide(y0(1:Fs*duration),Fs,fmin,fmax,500);
    y1 = fftFilterSingleSide(y1(1:Fs*duration),Fs,fmin,fmax,500);
    y2 = fftFilterSingleSide(y2(1:Fs*duration),Fs,fmin,fmax,500);
    y3 = fftFilterSingleSide(y3(1:Fs*duration),Fs,fmin,fmax,500);
    y4 = fftFilterSingleSide(y4(1:Fs*duration),Fs,fmin,fmax,500);
    y5 = fftFilterSingleSide(y5(1:Fs*duration),Fs,fmin,fmax,500);
    y6 = fftFilterSingleSide(y6(1:Fs*duration),Fs,fmin,fmax,500);
    y7 = fftFilterSingleSide(y7(1:Fs*duration),Fs,fmin,fmax,500);
    
    %inputSig = [y0(1:duration*Fs); y1(1:duration*Fs); y2(1:duration*Fs); y3(1:duration*Fs); y4(1:duration*Fs); y5(1:duration*Fs); y6(1:duration*Fs); y7(1:duration*Fs)]; % P by L (L is number of samples)
    rx = [y0(1:duration*Fs) y1(1:duration*Fs) y2(1:duration*Fs) y3(1:duration*Fs) y4(1:duration*Fs) y5(1:duration*Fs) y6(1:duration*Fs) y7(1:duration*Fs)]; % P by L (L is number of samples)
end

% beamformer1 = phased.MVDRBeamformer('SensorArray',array,...
%     'SampleRate',Fs,'PropagationSpeed',c,'OperatingFrequency',fc,...
%     'Direction',incidentAngle,'WeightsOutputPort',true);

% beamformer = phased.PhaseShiftBeamformer('SensorArray',array,...
%     'PropagationSpeed',c,'OperatingFrequency',fc,...
%     'Direction',incidentAngle,'WeightsOutputPort',true);


% beamformer = phased.TimeDelayBeamformer('SensorArray',array,...
%     'SampleRate',Fs,'PropagationSpeed',c,...
%     'Direction',incidentAngle,'WeightsOutputPort',true);

steervec = phased.SteeringVector('SensorArray',array,...
    'PropagationSpeed',c);
beamformer = phased.LCMVBeamformer('Constraint',steervec(fc,incidentAngle),'DesiredResponse',1,'WeightsOutputPort',true);

[y,w] = beamformer(rx);


figure;
plot(t,real(rx(:,3)),'r:',t,real(y),'b')
if (use_sim)
    xlim([0 0.001]);
end
xlabel('Time')
ylabel('Amplitude')
legend('Original','Beamformed')