clc
clear
tic
%% Reading
global Tracedata
S = dir('*.dat');
for k = 1:numel(S)
Headerdata=Seg2FileReader(S(k).name);
numberofgeophones=Headerdata.NumberOfTraces;
numberofsamples=Headerdata.TraceDescription.NumberOfSamples;

descaling_factor=zeros(numberofgeophones,1);
Geophonedata=zeros(numberofsamples,numberofgeophones);
Geophonelocation=zeros(numberofgeophones,1);

for i=1:numberofgeophones
var1=Headerdata.TraceDescription(i).TextData.DESCALING_FACTOR;
var2=split(var1,'');
descaling_factor(i)=str2double(var2(1));
dt=Headerdata.TraceDescription(i).TextData.SAMPLE_INTERVAL;
Geophonedata(:,i)=Tracedata(:,i)*descaling_factor(i);
Geophonelocation(i)=Headerdata.TraceDescription(i).TextData.RECEIVER_LOCATION;
end
Sourcelocation=Headerdata.TraceDescription(1).TextData.SOURCE_LOCATION;  

% if abs(Geophonelocation(1)-Sourcelocation)>abs(Geophonelocation(numberofgeophones)-Sourcelocation)
%     reverse=1;
%     Geophonedata=flip(Geophonedata,2);
% else
%     reverse=0;
% end

T=(dt:dt:(numberofsamples)*dt)';

%% Saving the data
filename=split(S(k).name,'.');
filename=cell2mat(strcat(filename(1),'.txt'));
writematrix([T Geophonedata], filename)
toc
end

% %% The below portion is only for monitoring and can be commented
% %% Display data
% Generaldata=Headerdata.FileHeader.TextData;
% Date=Generaldata.ACQUISITION_DATE;
% Time=Generaldata.ACQUISITION_TIME;
% Time=[Date '__' Time]
% Units=Generaldata.UNITS
% Numberofstacks=Headerdata.TraceDescription(1).TextData.STACK
% Sourcedistance=min(abs(Sourcelocation-Geophonelocation))
% Geophonespacing=abs(Geophonelocation(1)-Geophonelocation(2))
% Sampleinterval=dt
% 
% %% Plotting
% scalingfactor=5;
% Timeoffset=0;
% for i=1:numberofgeophones
%     plot(scalingfactor*Geophonedata(:,i)./max(max(abs(Geophonedata)))+(i)*Geophonespacing,T+Timeoffset)
%     hold on 
% end
% set(gca, 'YDir','reverse')
% ylim([0 max(T)])
% hold off
