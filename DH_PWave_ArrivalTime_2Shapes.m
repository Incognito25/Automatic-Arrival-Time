
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
    
ResultMatrix = [string('Depth (m)'), string('P-Wave File Name'), string('Window Width'), string('R2'), string('Arrival Time (ms)'), string('Channel Number')];

S = dir('*.txt');
for i = 1:numel(S)
    
%%

%     %Take the channel number as input 
%     prompt = {'Enter the Channel Number'};
%     dlgtitle = 'Channel Number';
%     dims = [1 50];
%     definput = {'1'};
%     opts.Interpreter = 'tex';
%     ChNo = inputdlg(prompt,dlgtitle,dims,definput,opts);
%     ChNo = str2num(cell2mat(ChNo));
    
    
 % Take the depth as input 
    prompt = {'Enter the Depth'};
    dlgtitle = 'Depth';
    dims = [1 50];
    definput = {'1'};
    opts.Interpreter = 'tex';
    Depth = inputdlg(prompt,dlgtitle,dims,definput,opts);
    Depth = str2num(cell2mat(Depth));
    
    
% Read the data from the Excel file
[num, ~, ~] = xlsread('DH E28 Data Sheet', 'Data Sheet');

% Extract the depth, noise, and P Wave columns
depth_values = num(:, 1);
p_wave_values = num(:, 2);

% Find the indices of the rows where the depth column matches the desired depth value
indices = find(depth_values == Depth);

% Extract the values of the noise and P Wave columns for those rows
PWaveFileName = string(p_wave_values(indices));


         PWaveFile = load(PWaveFileName +'.txt');

         time = PWaveFile(:,1)*1000; % 1000 is multiplied to convert seconds to miliseconds
         signal = -PWaveFile(:,ChNo+1)./max(abs(PWaveFile(:,ChNo+1)));

%% Make the wave start from zero

signal = signal - signal(1);
%%

% Invert the wave
figure
plot(time,signal);
title('Choose Valley or Peak')
ylabel('Amplitude')
xlabel('Time (ms)')
set(gcf, 'Position', get(0, 'Screensize'));
xlim([-5 100])

% Add textbox to the upper right hand side of the plot
annotation('textbox', [0.72 0.75 0.16 0.13], 'String', {'Instructions:', 'Draw a rectangle around any negative values to invert the wave, or vice versa.'}, 'EdgeColor', 'black','FontWeight','bold','FontSize',16);

% Draw the Peak Rectangle
xPeakRect = getrect;
xPeak1 =xPeakRect(1);
xPeak2 = xPeak1 + xPeakRect(3);
% Correcting xPeak1
diff = abs(time-xPeak1);
xPeakIndex1 = find((diff-min(diff))==0);
xDelay1 = time(xPeakIndex1);
% Correcting xPeak2
diff = abs(time-xPeak2);
xPeakIndex2 = find((diff-min(diff))==0);
xPeak2 = time(xPeakIndex2);

% Find the Delay from the drawn rectangle
timeVector = time(xPeakIndex1:xPeakIndex2);
signalVector = signal(xPeakIndex1:xPeakIndex2);
signs = sign(signalVector);

% Check if all of the signs are positive
if all(signs == 1)
    % The region is in the positive side of the waveform
    signal = signal;
elseif all(signs == -1)
    % The region is in the negative side of the waveform
    signal = -signal;
else
    % The region contains both positive and negative samples
end

pause(1)
%%
global x y uifig xdata
x = time;
y = signal;
data = [x y];
x1=1;    % x and y coordinates for the initial box
x2=6;
y1=-1;
y2=1;

dragdemo(x1,y1,x2,y2);


uiwait(uifig)
[~,ind1]=min(abs(x-min(xdata)));
[~,ind2]=min(abs(x-max(xdata)));
newdata=data(ind1:ind2,:);
figure
plot(newdata(:,1),newdata(:,2),'b');

% newTime = newdata(:,1);
% newSignal = newdata(:,2);

% Make the cropped peak start from zero
% newSignalRelative = newSignal-min(newSignal);

RelativeSignal = signal - min(signal);

%%
% Finding the Arrival Times with a Curve Fitting Technique

% Curve Fitting, Using Model 2: Two overlapping Gaussians (slower but fits wider range of peak shapes)
R_sqr_selected = 0;
cp_selected = 0;
ArrivalTime = 0;
WindowWidth = newdata(end,1)-newdata(1,1);
%To  create proper offset
if WindowWidth >1
    WindowWidth = floor(newdata(end,1)-newdata(1,1));
end

[~,max_idx] = max(newdata(:,2));
xPeak = newdata(max_idx,1);
xMid = newdata(round(end/2),1);

Offset = WindowWidth/2;
if Offset > 1
    Offset = 1;
end

if WindowWidth ==0
    WindowWidth = 1;
end

for cp = xPeak-Offset:0.05:xPeak+Offset;

%     Shape = 1; %Gaussian
% 
% [FitResults,GOF,baseline,coeff,residuals,xi,yi]=peakfit([time,RelativeSignal],cp,WindowWidth,3,Shape);
% 
% % Check if Shape is 1 and R2 is zero
% if (Shape == 1) && (GOF(2) <= 0)
%     Shape = 2; % Change Shape to 2
%     [FitResults,GOF,baseline,coeff,residuals,xi,yi] = peakfit([time,RelativeSignal],cp,WindowWidth,3,Shape); % Fit again with new Shape value
% end

    NoOfPeaks = 3; %No. Of Peaks to be fitted
        Shape = 1; %Gaussian

    while Shape <= 48

        try
            [FitResults,GOF,baseline,coeff,residuals,xi,yi]=peakfit([time,RelativeSignal],cp,WindowWidth,3,Shape);

            % Check if Shape is 1 and R2 is zero
            if (Shape == 1) && (GOF(2) == 0)
                Shape = 2; % Change Shape to 2
                continue; % Go back to the start of the loop
            end

            % Check if GOF(2) is less than or equal to 0 for current Shape
            if GOF(2) <= 0
                Shape = Shape + 1; % Change to next Shape
                continue; % Go back to the start of the loop
            end

            % If GOF(2) is greater than 0 for current Shape, break out of the loop
            break;

        catch
            % If there is an error in peakfit function, change to next Shape value
            Shape = Shape + 1;
        end
    end
    
    % Do something with the result
    if Shape <= 48
        disp(['Peak found at ', num2str(cp), ' with Shape = ', num2str(Shape)]);
    else
        disp(['No peak found at ', num2str(cp)]);
    end
        sizeresults=size(FitResults);
        if sizeresults(1)==1
            sy=yi; % total model for single-peak models
        else
            sy=sum(yi); % total model for multiple-peak models
        end
            maxsy=max(sy); % peak height maximum
            ymin=0.01.*maxsy; % cut-off amplitude
        % The next 3 statemens use theval2ind(x,val) function, which returns 
        % the index of the element of vector x that is closest to val.
        centerindex=val2ind(sy,maxsy); % point index number at peak center
        centerpoint=xi(val2ind(sy,sy(centerindex)));
        startpoint=xi(val2ind(sy(1:centerindex),ymin));
        R_sqr = GOF(2);
        if R_sqr > R_sqr_selected
            R_sqr_selected = R_sqr;
            cp_selected = cp;
            ArrivalTime = startpoint;
        else
            break
        end
end
R_sqr_selected;
ArrivalTime;
pause(1);
plotname = strcat('DH_',PWaveFileName,'_CURVE_FIT_',num2str(ChNo));
print(gcf, '-dtiff', '-r600', plotname)
%%
% Plot with the marked Arrival Time
figure('Name','Arrival Time','NumberTitle','on')
plot(time,signal,'LineWidth', 2.5, 'Color', 'blue');
hold on
xline(ArrivalTime,'--r','LineWidth', 2, 'Color', 'red');
hold off
title('Marked Arrival Time')
set(gca, 'FontWeight', 'bold')
% set(gcf, 'Position', get(0, 'Screensize'));
xlabel('Time (ms)')
ylabel('Amplitude')
xlim([-5 75])
ylim([min(signal)-0.1 max(signal)+0.1])
% dim = [0.15 0.8 0.1 0.1]; %left
dim = [0.57 0.85 0.05 0.05]; %right
str = {['DOWNHOLE SURVEY'],['-----------------------------------'],['File Name = ',num2str(PWaveFileName)],['Arrival Time = ',num2str(ArrivalTime),' ms'],['Window Width = ',num2str(WindowWidth)],['R^{2} = ',num2str(R_sqr_selected)],['Depth = ',num2str(Depth),' meters']};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
plotname = strcat('DH_',PWaveFileName,'_ARRIVAL_TIME_',num2str(ChNo));
print(gcf, '-dtiff', '-r600', plotname)
pause(2)

close all

 ResultMatrix = [ResultMatrix; [Depth PWaveFileName WindowWidth R_sqr  ArrivalTime ChNo]]
end
ResultMatrix;