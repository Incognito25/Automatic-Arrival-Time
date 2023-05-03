
clear
clc

% Take the channel number as input 
    prompt = {'Enter the Channel Number'};
    dlgtitle = 'Channel Number';
    dims = [1 50];
    definput = {'1'};
    opts.Interpreter = 'tex';
    ChNo = inputdlg(prompt,dlgtitle,dims,definput,opts);
    ChNo = str2num(cell2mat(ChNo));
    
ResultMatrix = [string('Depth of Source(m)'),string('Depth of Receiver(m)'), string('S+ Wave File Name'), string('S- Wave File Name'), string('Arrival Time (ms)')];

S = dir('*.txt');
for i = 1:numel(S)
    
  % Take the depth as input 
    prompt = {'Enter the Depth of Source'};
    dlgtitle = 'Depth';
    dims = [1 50];
    definput = {'1'};
    opts.Interpreter = 'tex';
    DepthOfSource = inputdlg(prompt,dlgtitle,dims,definput,opts);
    DepthOfSource = str2num(cell2mat(DepthOfSource));
    
    
% Read the data from the Excel file
[num, ~, ~] = xlsread('CH E14 Data Sheet', 'Data Sheet');

% Extract the depth, noise, and P Wave columns
depth_source_values = num(:, 1);
depth_receiver_values = num(:, 2);
SPlus_values = num(:, 4);
SMinus_values = num(:, 5);

% Find the indices of the rows where the depth column matches the desired depth value
indices = find(depth_source_values == DepthOfSource);


% Extract the values of the noise and P Wave columns for those rows
DepthOfSource = string(depth_source_values(indices));
DepthOfReceiver = string(depth_receiver_values(indices))
SPlusFileName = string(SPlus_values(indices));
SMinusFileName = string(SMinus_values(indices));

SPlusFileName = 'ST' + SPlusFileName;
SMinusFileName = 'ST' + SMinusFileName;

         yPlusFile = load(SPlusFileName +'.txt');
         xS = yPlusFile(:,1)*1000; % 1000 is multiplied to convert seconds to miliseconds
         yPlus = yPlusFile(:,ChNo+1)./max(abs(yPlusFile(:,ChNo+1)));
         
         
         yMinusFile = load(SMinusFileName +'.txt');
         yMinus = yMinusFile(:,ChNo+1)./max(abs(yMinusFile(:,ChNo+1)));
         
%%
minSize = min(numel(yPlus),numel(yMinus));
xS = xS(1:minSize);
yPlus = yPlus(1:minSize);
yMinus = yMinus(1:minSize);

%%
% Plot of the Original Data
figure ('Name','Original Data','NumberTitle','on')
plot(xS,yPlus,xS,yMinus)
set(gcf, 'Position', get(0, 'Screensize'));
legend('yPlus','yMinus');
title('Original Data');
xlim([0 100]);
plotname = strcat('ORIGINAL_DATA_',SPlusFileName,'_',SMinusFileName);
print(gcf, '-dtiff', '-r600', plotname);
%%
% Use a lowpass filter when the waves do not intersect at the required point (maybe 300Hz or so)
%          SamplingFrequency = 1000/(xS(2)-xS(1)); % divided by 1000 since time is in miliseconds
%          yPlusAfterSmooth = lowpass(yPlus,200,SamplingFrequency);
%          yMinusAfterSmooth = lowpass(yMinus,200,SamplingFrequency);
%          y = lowpass(x,fpass,fs) specifies that x has been sampled at a rate of fs hertz. fpass is the passband frequency of the filter in hertz
% The sampling rate of the data is 0.125 mili-seconds , thus 1/(0.125 * 10^-3) = 8000 Hz


% Using Smooth Function 
% yPlusSmooth = smooth(xS,yPlus,0.5,'rloess');
% yMinusSmooth = smooth(xS,yMinus,0.5,'rloess');
% 
% yPlusAfterSmooth = yPlus-yPlusSmooth;
% yMinusAfterSmooth = yMinus-yMinusSmooth;

% No Smooth
% yPlusAfterSmooth = yPlus;
% yMinusAfterSmooth = yMinus;

yPlusAfterSmooth = detrend(yPlus);
yMinusAfterSmooth = detrend(yMinus);

%%
yCommon = [];
xCommon = [];
%%
% Find the intersection of the curve

% Considering a 'Linear Bézier Curve' between two points
% Solve the system of linear equations with matrices
for i = 1:length(xS)-1
    A = [xS(i+1,:)-xS(i,:),0,-1,0;0,xS(i+1,:)-xS(i,:),-1,0;yPlusAfterSmooth(i+1,:)-yPlusAfterSmooth(i,:),0,0,-1;0,yMinusAfterSmooth(i+1,:)-yMinusAfterSmooth(i,:),0,-1];
    B = [-xS(i);-xS(i);-yPlusAfterSmooth(i);-yMinusAfterSmooth(i)];
% Solving for AT=B, thus T = inv(A)*B
% T = [t;u;xCommon;yCommon], where t & u are Bézier parameters
    T = A\B;              % T = inv(A)*B;
    if T(1,:)<=1 && T(1,:)>=0 && T(2,:)<=1 && T(2,:)>=0
        xCommon(i) = T(3,:);
        yCommon(i) = T(4,:);
    end
end

yCommon(yCommon==0) = NaN;
xCommon(xCommon==0) = NaN;

yCommon = yCommon';
xCommon = xCommon';

figure ('Name','Intersections','NumberTitle','on')
plot(xS,yPlusAfterSmooth,'LineWidth',2); 
hold on; 
plot(xS,yMinusAfterSmooth,'LineWidth',2); 
plot(xCommon,yCommon,'sr','LineWidth',2); 
title('Choose the Intersection Point');
set(gca, 'FontWeight', 'bold');
xlabel('Time (ms)');
ylabel('Amplitude');
xlim([-5 75]);
ylim([-1.1 1.1]);

plotname = strcat('INTERSECTIONS_',SPlusFileName,'_',SMinusFileName);
print(gcf, '-dtiff', '-r600', plotname);


%%
% Choosing the intersection point

% Draw the Common Rectangle
xCommonRect = getrect;
xCommon1 =xCommonRect(1);
xCommon2 = xCommon1 + xCommonRect(3);
% Correcting xCommon1
diff = abs(xS-xCommon1);
xCommonIndex1 = find((diff-min(diff))==0);
xCommon1 = xS(xCommonIndex1);
% Correcting xCommon2
diff = abs(xS-xCommon2);
xCommonIndex2 = find((diff-min(diff))==0);
xCommon2 = xS(xCommonIndex2);

% Find the Intersection from the drawn rectangle
xCommonIndex = find(xCommon > xCommon1 & xCommon < xCommon2);
xCommonSelected = xCommon(xCommonIndex);

ArrivalTime = xCommonSelected;
ArrivalTime = max(ArrivalTime);
%%
% Plot with the marked Arrival Time
figure('Name','Arrival Time','NumberTitle','on')
plot(xS,yPlus,'LineWidth',2);
hold on
plot(xS,yMinus,'LineWidth',2); 
xline(ArrivalTime,'--r','LineWidth',1.7); 
title('Marked Arrival Time');
set(gca, 'FontWeight', 'bold');
xlabel('Time (ms)');
ylabel('Amplitude');
xlim([-5 75]);
ylim([-1.1 1.1]);


dim = [0.57 0.85 0.05 0.05];
str = {['CROSSHOLE SURVEY'],['--------------------------------------'],['S-Plus File Name = ',num2str(SPlusFileName)],['S-Minus File Name = ',num2str(SMinusFileName)],['Arrival Time = ',num2str(ArrivalTime),' ms']};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',8);
plotname = strcat('CH_',SPlusFileName,'_',SMinusFileName,'_ARRIVAL_TIME');
print(gcf, '-dtiff', '-r600', plotname)
pause(1)

close all

ResultMatrix = [ResultMatrix; [DepthOfSource DepthOfReceiver SPlusFileName SMinusFileName ArrivalTime]]
end