
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
    
ResultMatrix = [string('Depth (m)'), string('S+ Wave File Name'), string('S- Wave File Name'), string('Arrival Time (ms)')];

S = dir('*.txt');
for i = 1:numel(S) %
 %%
 % Take the depth as input 
    prompt = {'Enter the Depth'};
    dlgtitle = 'Depth';
    dims = [1 50];
    definput = {'1'};
    opts.Interpreter = 'tex';
    Depth = inputdlg(prompt,dlgtitle,dims,definput,opts);
    Depth = str2num(cell2mat(Depth));
    
    
% Read the data from the Excel file
[num, ~, ~] = xlsread('DH_16th_December_2022.xlsx', 'Data_Sheet');

% Extract the depth, noise, and P Wave columns
depth_values = num(:, 1);
SPlus_values = num(:, 4);
SMinus_values = num(:, 5);

% Find the indices of the rows where the depth column matches the desired depth value
indices = find(depth_values == Depth);

% Extract the values of the noise and P Wave columns for those rows
SPlusFileName = string(SPlus_values(indices));
SMinusFileName = string(SMinus_values(indices));


         yPlusFile = load(SPlusFileName +'.txt');

         xS = yPlusFile(:,1)*1000; % 1000 is multiplied to convert seconds to miliseconds
         yPlus = yPlusFile(:,ChNo+1)./max(abs(yPlusFile(:,ChNo+1)));
         
         
         yMinusFile = load(SMinusFileName +'.txt');
         yMinus = yMinusFile(:,ChNo+1)./max(abs(yMinusFile(:,ChNo+1)));
%%
% Plot of the Original Data
figure ('Name','Original Data','NumberTitle','on')
plot(xS,yPlus,xS,yMinus)
set(gcf, 'Position', get(0, 'Screensize'));
legend('yPlus','yMinus');
title('Original Data');
xlim([0 100]);
pause(3)

%%

% Use a lowpass filter when the waves do not intersect at the required point (maybe 300Hz or so)
%          yPlus = lowpass(yPlus,250,8000);
%          yMinus = lowpass(yMinus,250,8000);
%          y = lowpass(x,fpass,fs) specifies that x has been sampled at a rate of fs hertz. fpass is the passband frequency of the filter in hertz
% The sampling rate of the data is 0.125 mili-seconds , thus 1/(0.125 * 10^-3) = 8000 Hz


% Using Smooth Function 
yPlusSmooth = smooth(xS,yPlus,0.1,'rloess');
yMinusSmooth = smooth(xS,yMinus,0.1,'rloess');

yPlusAfterSmooth = yPlus-yPlusSmooth;
yMinusAfterSmooth = yMinus-yMinusSmooth;

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

figure ('Name','Cropped Data','NumberTitle','on')
plot(xS,yPlusAfterSmooth,xS,yMinusAfterSmooth,xCommon,yCommon,'sr')
title('Choose the Intersection Point')
set(gcf, 'Position', get(0, 'Screensize'));
xlim([0 100])


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
plot(xS,yPlus,xS,yMinus,'LineWidth', 3);
hold on
xline(ArrivalTime,'--r','LineWidth', 2, 'Color', 'red');
hold off
title('Marked Arrival Time')
set(gcf, 'Position', get(0, 'Screensize'));
xlabel('Time (ms)')
ylabel('Amplitude')
xlim([-5 70])
dim = [0.15 0.5 0.3 0.3];
str = {['DOWNHOLE SURVEY'],['--------------------------------------'],['S-Plus File Name = ',num2str(SPlusFileName)],['S-Minus File Name = ',num2str(SMinusFileName)],['Arrival Time = ',num2str(ArrivalTime),' ms']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
plotname = strcat('DH_',SPlusFileName,'_',SMinusFileName,'_ARRIVAL_TIME');
saveas(gcf,plotname,'jpg')
pause(3)


ResultMatrix = [ResultMatrix; [Depth SPlusFileName SMinusFileName ArrivalTime]]
end
ResultMatrix;
close all
