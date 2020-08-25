function varargout = sineFitDemo(varargin)
%sineFitDemo
%Demo GUI for fitting noisy sine curves.
%Explanation
%   Input: all values are floating point, except 'num samples' must be integer!
%      'Number of periods' influences x end.
%         The number of periods may be less than 1, tested down to 0.1 periods.
%      'Samples per period' must be larger than 2.0. Interacts autom. with 'num samples'.
%      'Num samples' must be integer! Interacts autom. with 'Samples per period'
%          Minimum 4 samples in order to detect parameters.
%      'signal to noise': Noise is applied on the sine curve excluding offset
%         and is normaly in dB, depending on your Matlab version.
%   Output: Estimate of the noise free sine curve parameters.
%      MSE is taken from nlinfit
%         and is different to the MSE of the noise free function versus noisy input.
%   Best regression: The parameters of the clean sine are used as intial data for nlinfit.m.
%      If the result is different to the clean parameters, another sine fits the noisy sine better.
%   FFT: Parameters at the peak of the FFT. Offs. is the mean value of the input.
%   Process time: Time used by sineFit.m.
%   Run: Calculate the output for the input values.
%   Other noise: Run same clean sine with same SNR, but different random noise added.
%   Run xy.mat: runs your y(x), this must be in a 'xy.mat' file like:
%      x=[1,2,3,...]
%      y=[1,2,1,...]
%      Optional, if you know the noise free parameters add this variable:
%           paramsClean=[1,2,3,4];%[offset,amp,f,phase]
%           if you know in addition SNR:
%           paramsClean=[1,2,3,4,33];%[offset,amp,f,phase,SNR]
%           SNR is only for information.
%      see example file 'xy.mat'.
%    'Save In/Out': Save the input and output data to SineParamsInOut.mat.
%       x and y values of noisy input
%       paramsClean: Parameters of noise free sine:
%           [offset, amplitude, frequency, phase, SNR]
%       SineParamsOut: Estimated parameters of noisy sine:
%           [offset, amplitude, frequency, phase, MSE]
%       You may rename it to xy.mat and click on 'Run xy.mat'.
%       SineParamsOut is also present in the workspace.
%Definition of 10% error:
%	 The amplitude or the frequency is wrong,
%  if the deviation is more than 10% compared to the noise free sine.
%	 The offset can be zero and is considered to be wrong,
%  if it deviates by more than 10% of the expected offset
%  and more than 10% of the expected amplitude.
%	 The phase is wrong if it deviates by more than 10% of 2pi.
%Definition of 1% error:
%  Like above, but replace 10% with 1%.

% Last Modified by GUIDE v2.5 25-Oct-2019 09:55:37
%Author: Peter Seibold

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @sineFitDemo_OpeningFcn, ...
  'gui_OutputFcn',  @sineFitDemo_OutputFcn, ...
  'gui_LayoutFcn',  [] , ...
  'gui_Callback',   []);
if nargin && ischar(varargin{1})
  gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
  [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
  gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before sineFitDemo is made visible.
function sineFitDemo_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);% Update handles structure
CreateDrawSin(0,handles);
cmdRun_Callback(hObject, eventdata, handles);

% --- Outputs from this function are returned to the command line.
function varargout = sineFitDemo_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function editoffset_Callback(hObject, eventdata, handles)
CreateDrawSin(1,handles);

function editAmplitude_Callback(hObject, eventdata, handles)
CreateDrawSin(2,handles);

function editf_Callback(hObject, eventdata, handles)
CreateDrawSin(3,handles);

function editPhaseshift_Callback(hObject, eventdata, handles)
CreateDrawSin(4,handles);

function txtNumPeriods_Callback(hObject, eventdata, handles)
CreateDrawSin(5,handles);

function editsamplesperperiod_Callback(hObject, eventdata, handles)
CreateDrawSin(6,handles);

function editNumSamples_Callback(hObject, eventdata, handles)
CreateDrawSin(7,handles);

function editNoise_Callback(hObject, eventdata, handles)
CreateDrawSin(8,handles);

function edittstart_Callback(hObject, eventdata, handles)
CreateDrawSin(9,handles);

% --- Executes on button press in cmdRun.
function cmdRun_Callback(hObject, eventdata, handles)
global x y paramsClean SineParamsOut FFTparamsOut
FillEditBoxesInput(handles)
tic;
SineParamsOut=sineFit(x,y);
set(handles.txtProcessTime,'string',[sprintf('%0.0f',toc*1000) ' ms']);
%get MSE of output
modelfun = @(paramc,x) paramc(1) + paramc(2) * sin(2*pi* paramc(3)*x+paramc(4));
opts = statset('nlinfit');opts.MaxIter=0;
warning('off','all');
[~,~,~,~,MSE] = nlinfit(x,y,modelfun,SineParamsOut);
warning('on','all');
SineParamsOut(5)=MSE;
set(handles.txtMSEout,'string',['MSE: ' num2str(MSE,3)]);
%plot sine
x2=x(1):min(0.002*(x(end)-x(1)),(x(2)-x(1))/10):x(end);
y2=SineParamsOut(1)+SineParamsOut(2)*sin(2*pi*SineParamsOut(3)*x2+SineParamsOut(4));%result
axes(handles.plotFig);
pN(1)=plot(x,y,'Color',[0.8,0,0]);%input
hold on;
plot(x,y,'k.');%input dot
pN(2)=plot(x2,y2,'b-');%result
xlabel('Time [s]')
extra=(x2(end)-x2(1))*0.03;
xlim([x2(1)-extra x2(end)+extra]);
grid on;
if length(paramsClean)>3
  xstart=x(1);
  xend=x(end);
  extra=(xend-xstart)*0.02;
  xClean=xstart-extra:min(0.002*(xend-xstart),(x(2)-x(1))/10):xend+extra;
  yClean=paramsClean(1)+paramsClean(2)*sin(2*pi*paramsClean(3)*xClean+paramsClean(4));
  pN(3)=plot(xClean,yClean,'g-');%plot clean sine
  uistack(pN(2),'top');%put result line on top
  legend([pN(1),pN(3),pN(2)],'input','expected','output','Location','best');
  %best regression: use input values of GUI for regression
  initial_params=paramsClean(1:4);
  warning('off','all');
  opts = statset('nlinfit');opts.MaxIter=1000;
  [BETA,~,~,~,MSE] = nlinfit(x,y,modelfun,initial_params);
  warning('on','all');
  BETA(4)=rem(BETA(4),2*pi);
  if BETA(4)<0;BETA(4)=BETA(4)+2*pi; end;
  txtBestR=[num2str(BETA(1),4) '+' num2str(BETA(2),4) '*sin(2*pi*' num2str(BETA(3),4)...
    '*x+' num2str(BETA(4),4) ')  MSE: ' num2str(MSE,3)];
  set(handles.txtBestRegression,'string',txtBestR);
else
  legend([pN(1),pN(2)],'input','output','Location','best');
  set(handles.txtBestRegression,'string',' ');
end
hold off;
assignin('base', 'SineParamsOut', SineParamsOut);%you may delete this statement
PlotFFT(handles);
TxtBoxesErrorColor(handles);
%optional, x,y,paramsClean to workspace:
% assignin('base','paramsClean',paramsClean);
% assignin('base','xNoisy',x);
% assignin('base','yNoisy',y);

% --- Executes on button press in cmdOtherNoise.
function cmdOtherNoise_Callback(hObject, eventdata, handles)
global paramsClean
if length(paramsClean)<5
  set(handles.txtMessage,'String',...
    'ERROR: Not enogh parameters of clean sine!',...
    'Foregroundcolor',[.8 0 0]);
  return;
end
CreateDrawSin(0,handles);
cmdRun_Callback(hObject, eventdata, handles);

function cmdRunMyXY_Callback(hObject, eventdata, handles)
global paramsClean
TxtOutClear(handles);
set(handles.txtMessage,'String','');
%check if file exists
if exist(fullfile(pwd, 'xy.mat'), 'file')==2
  paramsClean=NaN;%overwritten if existing in 'xy.mat'
  load(fullfile(pwd, 'xy.mat'));
  cmdRun_Callback(hObject, eventdata, handles);
else
  set(handles.txtMessage,'String',...
    'ERROR: xy.mat not found in current directory!',...
    'Foregroundcolor',[.8 0 0]);
end

% --- Executes on button press in cmdSaveInOut.
function cmdSaveInOut_Callback(hObject, eventdata, handles)
global x y paramsClean SineParamsOut
save('SineParamsInOut.mat','x','y', 'paramsClean', 'SineParamsOut');

function cmdHelp_Callback(hObject, eventdata, handles)
help sineFitDemo;

function CreateDrawSin(Flag,handles)
%runs if input values (offset, amplitude, etc.) changed by user or with cmdRun
global x y paramsClean
TxtOutClear(handles);%Clear output numbers
set(handles.txtMessage,'String','');%Clear message box
cla(handles.plotFigFFT,'reset');%Clear FFT plot
offset=str2double(get(handles.editoffset,'string'));
amplitude=str2double(get(handles.editAmplitude,'string'));
f=str2double(get(handles.editf,'string'));
Phi0=str2double(get(handles.editPhaseshift,'string'));
numPeriods=str2double(get(handles.txtNumPeriods,'string'));
samplesPperiod=str2double(get(handles.editsamplesperperiod,'string'));
numSamples=str2double(get(handles.editNumSamples,'string'));
NoiseLevel=str2double(get(handles.editNoise,'string'));
xstart=str2double(get(handles.edittstart,'string'));
%adjust samplesPperiod, NumSamples accordingly:
if Flag==5  && sum(isnan([numPeriods,samplesPperiod]))==0
  %Num periods changed by user
  %Keep samples per period if num periods >1, else keep num samples
  if numPeriods>=1
    numSamples=max([round(numPeriods*samplesPperiod),...
      ceil(2*numPeriods+eps(2*numPeriods)),4]);
  end
elseif Flag==6 && sum(isnan([numPeriods,samplesPperiod]))==0
  %Samples per period changed by user
  numSamples=max(round(numPeriods*samplesPperiod),3);
elseif Flag==7 && ~isnan(numSamples)
  %Num samples changed by user
  %numSamples must be an integer number!
  numSamples=max(round(numSamples),3);
end
if numSamples<=2*numPeriods || numSamples<4
  Ns=max(ceil(2*numPeriods+eps(2*numPeriods)),4);
  set(handles.txtMessage,'String',...
    ['WARNING: Not enough samples to detect parameters, num samples >= ' num2str(Ns)],...
    'Foregroundcolor',[.7 .4 0]);
end
samplesPperiod=numSamples/numPeriods;
set(handles.editsamplesperperiod,'String',num2str(samplesPperiod,6));
set(handles.editNumSamples,'String',num2str(numSamples));
%check if inputs are numbers:
if sum(isnan([offset,amplitude,f,Phi0,numPeriods,samplesPperiod,numSamples,NoiseLevel,xstart]))>0
  cla(handles.plotFig,'reset');%clear sine plot
  set(handles.txtMessage,'String',...
    'ERROR INPUT: INVALID NUMBERS PRESENT1!',...
    'Foregroundcolor',[.8 0 0]);
  return
end
%create sine acc. to user input
xend=xstart+numPeriods/f;
x=xstart:1/(f*samplesPperiod):xend;
x=x(1:end-1);
xend=x(end);
paramsClean=[offset,amplitude,f,Phi0,NoiseLevel];
set(handles.edittend,'string',num2str(xend,4));
set(handles.editNumSamples,'string',num2str(length(x),4));
ysin=amplitude*sin(2*pi*f*x+Phi0);
if NoiseLevel<99 % no noise if SNR>=99
  if exist('awgn.m','file')==2
    %communication_toolbox necessary
    ysin=awgn(ysin,NoiseLevel,'measured');
  elseif exist('randn.m','file')==2
    %Signal Processing Toolbox necessary
    ysin=ysin+amplitude*randn(size(ysin))/(10^(NoiseLevel*20));
  else
    ysin=ysin+amplitude*(rand(size(ysin))-0.5)*2/(10^(NoiseLevel*20));
  end
end
y=offset+ysin;
%plot created sine
extra=(xend-xstart)*0.02;% for longer green line
xClean=xstart-extra:min(0.002*(xend-xstart),(x(2)-x(1))/10):xend+extra;
yClean=offset+amplitude*sin(2*pi*f*xClean+Phi0);
axes(handles.plotFig);
pN=plot(x,y,'k.',xClean,yClean,'g-');%dot input, clean sine
hold on;
pN(3)=plot(x,y,'Color',[0.8,0,0]);% line input
legend([pN(3),pN(2)],'input','expected','Location','best');
hold off;
xlabel('Time [s]')
extra=(xend-xstart)*0.01;
xlim([xClean(1)-extra xClean(end)+extra]);
grid on;

function PlotFFT(handles)
global x y paramsClean SineParamsOut FFTparamsOut
numSamples=length(x);
Ts=x(2)-x(1);%sampling Period
offs=mean(y);%DC value
y_m=y-offs;
n = 128*2^nextpow2(numSamples);%heavy zero padding
Y = fft(y_m,n);%Y(f)
n2=n/2;
P2 = abs(Y/numSamples);
P1 = P2(1:n2+1);
P1(2:end-1) = 2*P1(2:end-1);
fs = (0:n2)/n/Ts;% f scale
[maxFFT,maxFFTindx]=max(P1);
fPeak=fs(maxFFTindx);% f at peak
%Phase
Phip=angle(Y(maxFFTindx))+pi/2;%sin(90°+alpha)=cos(betta), alpha=-betta
%phase not correct due to rectangular window effect. OK for more than 1 period
Phip=Phip-x(1)*fPeak*2*pi;%shift for phi at x=0
Phip=rem(Phip,2*pi);
if Phip<0
  Phip=Phip+2*pi;
end
% %Better estimate for offset:
omega=2*pi*fPeak;
offs=offs-maxFFT*(cos(omega*x(1)+Phip)-cos(omega*x(end)+Phip))/(omega*(x(end)-x(1)));
FFTparamsOut=[offs maxFFT fPeak Phip];
%plot FFT
axes(handles.plotFigFFT);
pFFTin=plot(fs,P1,'r-');
xlabel('Frequency [Hz]');
ylabel('Amplitude')
hold on;
pFFTmax=plot(fs(maxFFTindx),maxFFT,'r+','MarkerSize',12);
pFFTresult=plot(SineParamsOut(3),SineParamsOut(2),'b+','LineWidth',2,'MarkerSize',8);
plot([SineParamsOut(3),SineParamsOut(3)],[0,max(max(P1)*1.01,SineParamsOut(2))],'b-');
if length(paramsClean)>3
  %paramsClean present
  pFFTexp=plot(paramsClean(3),paramsClean(2),'gx','LineWidth',2);
  plot([paramsClean(3),paramsClean(3)],[0,max(maxFFT*1.01,paramsClean(2)*1.01)],'g-');
  hLeg=legend([pFFTin,pFFTexp,pFFTresult,pFFTmax],'Input',...
    ['Expected: ' num2str(paramsClean(2),3) ' | ' num2str(paramsClean(3),3) ' Hz'],...
    ['Result:     ' num2str(SineParamsOut(2),3) ' | ' num2str(SineParamsOut(3),3) ' Hz'],...
    ['max FFT:  ' num2str(maxFFT,3) ' | ' num2str(fs(maxFFTindx),3) ' Hz'],...
    'Location','best');
else
  hLeg=legend([pFFTin,pFFTresult,pFFTmax],'Input',...
    ['Result:       ' num2str(SineParamsOut(2),3) ' | ' num2str(SineParamsOut(3),3) ' Hz'],...
    ['max FFT:  ' num2str(maxFFT,3) ' | ' num2str(fs(maxFFTindx),3) ' Hz'],...
    'Location','best');
end
title(hLeg,'        amplitude | frequency','FontSize',8);
hold off;
grid on;

function FillEditBoxesInput(handles)
global x  paramsClean
if length(paramsClean)>3
  set(handles.editoffset,'string',num2str(paramsClean(1),4));
  set(handles.editAmplitude,'string',num2str(paramsClean(2),4));
  set(handles.editf,'string',num2str(paramsClean(3),4));
  set(handles.editPhaseshift,'string',num2str(paramsClean(4),4));
  numPeriods=(x(end)-x(1)+x(2)-x(1))*paramsClean(3);
  samplesperperiod=(length(x))/numPeriods;
  set(handles.txtNumPeriods,'string',num2str(numPeriods,6));
  set(handles.editsamplesperperiod,'string',num2str(samplesperperiod,6));
  if length(x)<=2*numPeriods || length(x)<4
    Ns=max(ceil(2*numPeriods+eps(2*numPeriods)),4);
    set(handles.txtMessage,'String',...
      ['WARNING: Not enough samples to detect parameters, num samples >= ' num2str(Ns)],...
      'Foregroundcolor',[.7 .4 0]);
  end
else
  set(handles.editoffset,'string','');
  set(handles.editAmplitude,'string','');
  set(handles.editf,'string','');
  set(handles.editPhaseshift,'string','');
  set(handles.txtNumPeriods,'string','');
  set(handles.editsamplesperperiod,'string','');
end
set(handles.editNumSamples,'string',num2str(length(x)));
set(handles.edittstart,'string',num2str(x(1),4));
set(handles.edittend,'string',num2str(x(end),4));
if length(paramsClean)==5
  set(handles.editNoise,'string',num2str(paramsClean(5),3));
else
  set(handles.editNoise,'string','');
end

function TxtOutClear(handles)
set(handles.txtOffset,'BackgroundColor',[1 1 1]);
set(handles.txtAmplitude,'BackgroundColor',[1 1 1]);
set(handles.txtf,'BackgroundColor',[1 1 1]);
set(handles.txtPhi0,'BackgroundColor',[1 1 1]);
set(handles.txtOffset,'string',' ');
set(handles.txtAmplitude,'string',' ');
set(handles.txtf,'string',' ');
set(handles.txtPhi0,'string',' ');
set(handles.txtBestRegression,'string',' ');
set(handles.txtFFToffs,'string',' ');
set(handles.txtFFTamp,'string',' ');
set(handles.txtFFTf,'string',' ');
set(handles.txtFFTphi,'string',' ');

function TxtBoxesErrorColor(handles)
global paramsClean SineParamsOut FFTparamsOut
if length(paramsClean)>3
  ErrorL1 = 0.01; % 1% error
  ErrorL10 = 0.1; % 10% error
  errorPi1 =2*pi*0.01; % 1% phase error
  errorPi10 = errorPi1*10; % 10% phase error
  green=[.9 1 .9];
  orange=[1 .95 .8];
  red=[1 .8 .8];
  OffClean=paramsClean(1);
  AmpClean=paramsClean(2);
  fClean=paramsClean(3);
  PhiClean=paramsClean(4);
  OffOut=SineParamsOut(1);
  AmpOut=SineParamsOut(2);
  fOut=SineParamsOut(3);
  PhiOut=SineParamsOut(4);
  %Offset
  set(handles.txtOffset,'string',num2str(OffOut(1),4));
  if abs(OffOut - OffClean) < ErrorL1 * AmpClean || abs(OffOut - OffClean) < ErrorL1 * abs(OffClean)
    set(handles.txtOffset,'BackgroundColor',green);
  elseif abs(OffOut - OffClean) < ErrorL10 * AmpClean || abs(OffOut - OffClean) < ErrorL10 * abs(OffClean)
    set(handles.txtOffset,'BackgroundColor',orange);
  else
    set(handles.txtOffset,'BackgroundColor',red);
  end
  %Amplitude
  set(handles.txtAmplitude,'string',num2str(AmpOut,4));
  if abs(AmpOut / AmpClean - 1) < ErrorL1
    set(handles.txtAmplitude,'BackgroundColor',green);
  elseif  abs(AmpOut / AmpClean - 1) < ErrorL10
    set(handles.txtAmplitude,'BackgroundColor',orange);
  else
    set(handles.txtAmplitude,'BackgroundColor',red);
  end
  %Frequency
  set(handles.txtf,'string',num2str(fOut,4));
  if abs(fOut / fClean - 1) < ErrorL1
    set(handles.txtf,'BackgroundColor',green);
  elseif abs(fOut / fClean - 1) < ErrorL10
    set(handles.txtf,'BackgroundColor',orange);
  else
    set(handles.txtf,'BackgroundColor',red);
  end
  %Phase
  set(handles.txtPhi0,'string',num2str(PhiOut,3));
  if PhiOut < 1 && PhiClean > 5
    PhiOut = PhiOut + 2 * pi;
  elseif PhiOut > 5 && PhiClean < 1
    PhiOut = PhiOut - 2 * pi;
  end
  if abs(PhiOut - PhiClean) < errorPi1
    set(handles.txtPhi0,'BackgroundColor',green);
  elseif abs(PhiOut - PhiClean) < errorPi10
    set(handles.txtPhi0,'BackgroundColor',orange);
  else
    set(handles.txtPhi0,'BackgroundColor',red);
  end
  % for FFT
  OffOut=FFTparamsOut(1);
  AmpOut=FFTparamsOut(2);
  fOut=FFTparamsOut(3);
  PhiOut=FFTparamsOut(4);
  %Offset
  set(handles.txtFFToffs,'string',num2str(OffOut(1),4));
  if abs(OffOut - OffClean) < ErrorL1 * AmpClean || abs(OffOut - OffClean) < ErrorL1 * abs(OffClean)
    set(handles.txtFFToffs,'BackgroundColor',green);
  elseif abs(OffOut - OffClean) < ErrorL10 * AmpClean || abs(OffOut - OffClean) < ErrorL10 * abs(OffClean)
    set(handles.txtFFToffs,'BackgroundColor',orange);
  else
    set(handles.txtFFToffs,'BackgroundColor',red);
  end
  %Amplitude
  set(handles.txtFFTamp,'string',num2str(AmpOut,4));
  if abs(AmpOut / AmpClean - 1) < ErrorL1
    set(handles.txtFFTamp,'BackgroundColor',green);
  elseif  abs(AmpOut / AmpClean - 1) < ErrorL10
    set(handles.txtFFTamp,'BackgroundColor',orange);
  else
    set(handles.txtFFTamp,'BackgroundColor',red);
  end
  %Frequency
  set(handles.txtFFTf,'string',num2str(fOut,4));
  if abs(fOut / fClean - 1) < ErrorL1
    set(handles.txtFFTf,'BackgroundColor',green);
  elseif abs(fOut / fClean - 1) < ErrorL10
    set(handles.txtFFTf,'BackgroundColor',orange);
  else
    set(handles.txtFFTf,'BackgroundColor',red);
  end
  %Phase
  set(handles.txtFFTphi,'string',num2str(PhiOut,3));
  if PhiOut < 1 && PhiClean > 5
    PhiOut = PhiOut + 2 * pi;
  elseif PhiOut > 5 && PhiClean < 1
    PhiOut = PhiOut - 2 * pi;
  end
  if abs(PhiOut - PhiClean) < errorPi1
    set(handles.txtFFTphi,'BackgroundColor',green);
  elseif abs(PhiOut - PhiClean) < errorPi10
    set(handles.txtFFTphi,'BackgroundColor',orange);
  else
    set(handles.txtFFTphi,'BackgroundColor',red);
  end
end

% --- Executes when SineFitDemoFig is resized.
function SineFitDemoFig_SizeChangedFcn(hObject, eventdata, handles)
FigPos=get(handles.SineFitDemoFig,'Position');
set(handles.uipanel,'Position',[1 FigPos(4)-105 FigPos(3)-1 105]);
