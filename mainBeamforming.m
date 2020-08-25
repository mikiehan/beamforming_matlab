radius=0.0383; % Radius (m) % 38.3 mm for matrix voice

%M # of theta scanning angles, N # of phi scanning angless, P number of mics
thetaScanningAngles = 0:(2*pi/360):2*pi;
phiScanningAngles = 90; 

M = numel(thetaScanningAngles);
N = numel(phiScanningAngles);
P=8; % Number of microphones
s=1; % Number of source signals
noise_var=0;


ebi = zeros(M, N, P); % steering vector;

xPos = [0, -38.13, -20.98, 11.97, 35.91, 32.81, 5.00, -26.57] ./ 38.3;
yPos = [0, 3.58, 32.04, 36.38, 13.32, -19.77, -37.97, -27.58] ./ 38.3;
zPos = zeros(1, P);
[gamma, l] = cart2pol(xPos, yPos);

folderName = '90';

[y0,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_0.wav'));
[y1,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_1.wav'));
[y2,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_2.wav'));
[y3,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_3.wav'));
[y4,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_4.wav'));
[y5,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_5.wav'));
[y6,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_6.wav'));
[y7,Fs] = audioread(strcat('/Users/profhan/Downloads/' , folderName , '/channel_7.wav'));

inputSig = [y0 y1 y2 y3 y4 y5 y6 y7];
inputSig = inputSig'; % P by L (L is number of samples)

c = 3.40e2; % Speed of sound (m/s)
fi = 19e3; % center frequency of FMCW from 18kHz - 20 kHz

%steeringVec_(j,0) = cos(2*PI*16000*((0.054*(cos((i*PI)/180.0) - cos( j*((360.0*PI)/(8.0*180.0)) - ((i*PI)/180.0))))/330.0));
 
% complex ebi
ebi = steeringVector(xPos, yPos, zPos, gamma, radius, true , fi, c, thetaScanningAngles, phiScanningAngles);

% real ebi 
%ebi = cos(2*pi*16000*((radius*(cos((angle*pi)/180) - cos(channel * ((360*pi)/(8 * 180)) - ((angle*pi)/180))))/330.0));
wi = weightingVectorMVDR(inputSig, ebi);
beamformed_Signal = zeros(M, length(y0));
for angle = 1:M
    weights = zeros(1, P);
    for channel = 1:P
        weights(channel) = wi(angle,1, channel);
    end
    beamformed_Signal(angle, :) = weights *  inputSig; 
end

%AoA index is 91
figure;
plot(abs(beamformed_Signal(91, :)));

%steering vector complex? vs real 
