% delayed estimation task with a single orientation stumulus in the center

function run_DE
Screen('Preference', 'SkipSyncTests', 1)

% shuffle the random seed
rng('shuffle');

% open window
windowPtr = Screen('OpenWindow',0);

% Screen information
screen.width = 40;
screen.number = max(Screen('Screens'));
[w, h]=Screen('WindowSize', screen.number);  % Screen resolution
screen.resolution = [w h];
screen.distance = 60;                      % distance between observer and Screen (in cm)
screen.angle = 2*(180/pi)*(atan((screen.width/2) / screen.distance)) ; % total visual angle of Screen
screen.ppd = screen.resolution(1) / screen.angle;  % pixels per degree
screen.bgcolor = [128 128 128];                % Screen bg color (grey)

% hide the cursor
HideCursor;

% fixation
x = w/2;                % fixation point (x)
y = h/2;                % fixation point (y)
fixcolor = [0 0 0];     % color of the fixation cross; default
armsize = 4;            % fixation cross size (each arm) 

% stimuli
settings.minaxis = .41 * screen.ppd;    % minor axis of the ellipse
settings.majaxis = .94 * screen.ppd;    % major axis of the ellipse
settings.tot_trials = 10;               % total # of trials
settings.stimtime = 100;                % stimulus time in ms
settings.delay = 1;                     % delay time btw two displays in second

% initialize
nextFlipTime = 0;                       % flip time for 'Flip' command
trialcnt = 0;                           % trial count
datamatrix = zeros(settings.tot_trials,19);      % data matrix

%% Screen 0: inter-trial display

Screen('FillRect',windowPtr,screen.bgcolor);
drawfixation(windowPtr, x, y, fixcolor, armsize);
Screen('Flip',windowPtr,nextFlipTime);
nextFlipTime = GetSecs + 1;

%% Screen 1: 1st display

Screen('FillRect',windowPtr,screen.bgcolor);
drawfixation(windowPtr, x, y, fixcolor, armsize);
ori = -90 + rand*180;           % orientation information of stimuli 

% draw stimuli
im = drawEllipse(settings.minaxis, settings.majaxis , ori, 200, 128);
texture = Screen('MakeTexture',windowPtr,im);
rect = CenterRectOnPoint([0 0 size(im')], x, y);
Screen('DrawTexture', windowPtr, texture ,[0 0 size(im')], rect);
Screen('Flip',windowPtr,nextFlipTime);

tic;
nextFlipTime = GetSecs + (settings.stimtime -12)/1000;

%% Screen 2: blank

Screen('FillRect',windowPtr,screen.bgcolor);
% drawfixation(windowPtr, x, y, fixcolor, armsize);
Screen('Flip',windowPtr,nextFlipTime);
STIMTIME1 = round(toc*1000);
nextFlipTime = GetSecs + settings.delay;

%% Screen 3: response

% 1. show ellipse with a random orientation
ori_probe = -90 + rand*180;            
im = drawEllipse(settings.minaxis, settings.majaxis , ori_probe, 200, 128);
texture = Screen('MakeTexture',windowPtr,im);
rect = CenterRectOnPoint([0 0 size(im')], x, y);
Screen('DrawTexture', windowPtr, texture ,[0 0 size(im')], rect);
Screen('Flip',windowPtr,nextFlipTime);

% 2. read mouse input
SetMouse(x, y);     % initialize mouse position

% initialize
respIdx = -1;
aborted = 0;
while (respIdx == -1) && ~aborted
% 3. change orientation according to the mouse
    [mouseX,mouseY,~] = GetMouse;                            % get mouse position
    if x==mouseX
        angle_value = 90;
    else
        angle_value = atan((mouseY-y)/(x-mouseX)) * 180/pi + 90; % get angle [0,180]
    end
    true_angle_value = angle_value;
    
    im = drawEllipse(settings.minaxis, settings.majaxis , angle_value, 200, 128);
    texture = Screen('MakeTexture',windowPtr,im);
    rect = CenterRectOnPoint([0 0 size(im')], x, y);
    Screen('DrawTexture', windowPtr, texture ,[0 0 size(im')], rect);
    Screen('Flip',windowPtr,nextFlipTime);
    
    % 4. press space bar to report#
    [~,~,keyCode] = KbCheck;
    if keyCode(KbName('space'))
%         Screen('DrawText',windowPtr,num2str(true_angle_value,3),x,y-100,0);
        respIdx=0;
    end
end

%% Screen 4: feedback
Screen('FillRect',windowPtr,screen.bgcolor);

% stimulus
im = drawEllipse(settings.minaxis, settings.majaxis , ori, 200, 128);
texture = Screen('MakeTexture',windowPtr,im);
rect = CenterRectOnPoint([0 0 size(im')], x, y-100);
Screen('DrawTexture', windowPtr, texture ,[0 0 size(im')], rect);
Screen('DrawText',windowPtr,'stimulus',x-150,y-120,0);

% response
im = drawEllipse(settings.minaxis, settings.majaxis , angle_value, 200, 128);
texture = Screen('MakeTexture',windowPtr,im);
rect = CenterRectOnPoint([0 0 size(im')], x, y+100);
Screen('DrawTexture', windowPtr, texture ,[0 0 size(im')], rect);
Screen('DrawText',windowPtr,'response',x-150,y+80,0);

nextFlipTime = GetSecs + .8;
Screen('Flip',windowPtr,nextFlipTime);

function im = drawEllipse(d1,d2,rot,fgcol,bgcol)

rot = -rot-90;  

d1 = round(d1);
d2 = round(d2);

rot = -rot/180*pi;

% make sure that d1 is the minor axis
if (d1>d2)
    d3=d1;
    d1=d2;
    d2=d3;
end

% draw ellipse
im = ones(2*d2,2*d2)*bgcol;
minX = -d2;
maxX = minX + 2*d2 - 1; 
[X Y] = meshgrid(minX:maxX,minX:maxX);
X_new = X * cos(rot) - Y * sin(rot);
Y = X * sin(rot) + Y * cos(rot);
X = X_new;

idx = (X.^2/(d1/2)^2 + Y.^2/(d2/2)^2)<1;
idx_low = (X.^2/((1.05*d1)/2)^2 + Y.^2/((1.025*d2)/2)^2)<1;
idx_super_low = (X.^2/((1.075*d1)/2)^2 + Y.^2/((1.05*d2)/2)^2)<1;

im(idx_super_low) = mean([fgcol bgcol bgcol]);
im(idx_low) = mean([fgcol bgcol]);
im(idx) = fgcol;

% rotate
% % im=imrotate(im,-rot);
% im(im==0)=bgcol;

% crop
while im(:,1)==bgcol
    im = im(:,2:end);
end
while im(1,:)==bgcol
    im = im(2:end,:);
end
while im(end,:)==bgcol
    im = im(1:end-1,:);
end
while im(:,end)==bgcol
    im = im(:,1:end-1);
end

function drawfixation(windowPtr, x, y, color, armsize)
Screen('DrawLine', windowPtr, color, x-armsize, y, x+armsize, y, 2);
Screen('DrawLine', windowPtr, color, x, y-armsize, x, y+armsize, 2);


