% function experiment(subjid,delaytime)
%
% Delay time should be specified in milliseconds
%
% E.g. experiment('RB',1000)  

function experiment(subjid,delaytime)

s = RandStream.create('mt19937ar','seed',sum(100*clock));
rng('shuffle')

if ~exist('subjid','var')
    clc;
    subjid = input('Your initials: ','s');
end

%-%-%-%-%-
%- INIT %-
%-%-%-%-%-
screen_width = 40;    % in cm (Dell@T115A: ~48cm; Dell@T101C: ~40 cm)
outputfile = ['output/' upper(subjid) '_' datestr(clock,30) '_' num2str(delaytime) '.mat'];

settings.maxdelay = 4000;  % maximum delay period (in ms)
% settings.intertrialtime = settings.maxdelay+500-settings.delay; % time between two trials
settings.intertrialtime = settings.maxdelay+500-delaytime; % time between two trials
settings.bglum = 30;
settings.gablambda = .25;
settings.gabsigma = .2;
settings.gabphase = 0;
settings.gabpeak = 0.8;   % between 0 and 1 (1=maximum contrast)
settings.bgdac = 128;
settings.stimtime = 90;
settings.delaytime = delaytime;
settings.stimecc = 5;
settings.setsizes = [2 4 6 8];
settings.nrep = 100;   % total number of repetitions of each set size (total #trials = length(setsizes)*nrep)
settings.breakTime = 30;

if delaytime>settings.maxdelay
    error('Chosen delay time is larger than settings.maxdelay!');
end

nTrials = length(settings.setsizes)*settings.nrep;
breakpoints = floor(nTrials*[.2 .4 .6 .8]);  % insert breaks after completeing 20%, 40%, 60%, and 80% of the trials

settings.kappa_prior = 0;

% screen info
screenNumber=0;
[w h]=Screen('WindowSize', screenNumber);  % screen resolution
screen_resolution = [w h];                 % screen resolution
screen_distance = 60;                      % distance between observer and screen (in cm)
screen_angle = 2*(180/pi)*(atan((screen_width/2) / screen_distance)) ; % total visual angle of screen
screen_ppd = screen_resolution(1) / screen_angle;  % pixels per degree
screen_fixposxy = screen_resolution .* [.5 .5]; % fixation position

% open screen
HideCursor;
gray=GrayIndex(screenNumber);
windowPtr = screen('OpenWindow',screenNumber,gray,[],32,2);

% show start screen
xpos = 100;
ypos = 20;
dy = 37;
screen('fillRect',windowPtr,128);
screen('TextSize',windowPtr,25);
screen('DrawText',windowPtr,['This experiment consists of ' num2str(nTrials) ' trials.'],xpos,ypos,[255 255 255]); ypos = ypos+2*dy;
screen('DrawText',windowPtr,'On each trial, a set of stimuli will be shown.',xpos,ypos,[255 255 255]); ypos = ypos+dy;
screen('DrawText',windowPtr,'After a brief delay, you will be asked to',xpos,ypos,[255 255 255]); ypos = ypos+dy;
screen('DrawText',windowPtr,'estimate the orienation of one of them. ',xpos,ypos,[255 255 255]); ypos = ypos+2*dy;
screen('DrawText',windowPtr,'The stimulus orientations are completely random.',xpos,ypos,[255 255 255]); ypos = ypos+dy;
screen('DrawText',windowPtr,'Here is an example set:',xpos,ypos,[255 255 255]); ypos = ypos+3*dy;
curry = ypos;
dx=60;
dy=60;
for ii=1:6
    currx=45;
    for jj=1:16  
        currx = currx+dx;
        randort = circ_vmrnd(0,settings.kappa_prior);
        randort = randort/pi*90;
        im = generate_gabor(randort,settings.gablambda*screen_ppd,settings.gabsigma*screen_ppd,settings.gabphase);
        im = ((settings.gabpeak*im*256)+256)/2;
        patchsize = size(im);
        screen('PutImage',windowPtr,im,centerRectOnPoint([0 0 size(im)],currx,curry));
        stimtex=screen('MakeTexture', windowPtr, im); 
    end
    curry=curry+dy;   
end
dy=37;
ypos = curry+dy;
screen('DrawText',windowPtr,'Good luck :-)',xpos,ypos,[255 255 255]); ypos = ypos+dy;
screen('DrawText',windowPtr,'Press any key to start',xpos,ypos,[60 200 60]);
screen('flip',windowPtr);
% im=screen('GetImage',windowPtr,[50 50 currx+30 curry+20]);
% imwrite(im,['example_' num2str(settings.kappa_prior) '.png'],'png');
waitForKey;
screen('FillRect', windowPtr, gray);
drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
screen('Flip', windowPtr);
screen('TextSize', windowPtr, 15);
waitSecs(1.5);

%-%-%-%-%-%-%-%-%-%
% Generate trials %
%-%-%-%-%-%-%-%-%-%
data.N = repmat(settings.setsizes,1,settings.nrep);  % 5 set sizes, 450 trials per session
data.N = data.N(randperm(length(data.N)));
for ii=1:length(data.N)
    data.stimvec{ii} = circ_vmrnd(zeros(1,data.N(ii)),settings.kappa_prior)/pi*90;
    data.targetidx(ii) = randi(data.N(ii));
    data.targetval(ii) = data.stimvec{ii}(data.targetidx(ii));
    data.startpos(ii) = randi(8);
end

% compute stimulus positions
angle=0;
for ii=1:8
    [x y] = pol2cart(angle,settings.stimecc*screen_ppd);
    posx(ii) = x+screen_fixposxy(1);
    posy(ii) = y+screen_fixposxy(2);
    angle = angle+2*pi/8;
end
 
%-%-%-%-%-%-%-%-%-%-%-%-%
%- LOOP THROUGH TRIALS %-
%-%-%-%-%-%-%-%-%-%-%-%-%
nextFlipTime = 0; % just to initialize...
aborted = 0;

for trialnr = 1:length(data.N)

    % create stimulus patches
    clear stimtex;
    for ii=1:data.N(trialnr)
        im = generate_gabor(data.stimvec{trialnr}(ii),settings.gablambda*screen_ppd,settings.gabsigma*screen_ppd,settings.gabphase);
        im = ((settings.gabpeak*im*256)+256)/2;
        patchsize(ii,:) = size(im);
        stimtex(ii)=screen('MakeTexture', windowPtr, im);
    end
    
    % SCREEN 1: FIXATION
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
    screen('flip',windowPtr,nextFlipTime);
    nextFlipTime = getSecs + .5;
    
    % SCREEN 2: STIMULUS
    screen('fillRect',windowPtr,settings.bgdac);    
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
    posidx = data.startpos(trialnr);
    for ii=1:data.N(trialnr)
        srcrect = [0 0 patchsize(ii,:)];        
        destrect = centerRectOnPoint(srcrect,posx(posidx),posy(posidx));
        screen('drawtexture',windowPtr,stimtex(ii),srcrect,destrect);
        if ii==data.targetidx(trialnr)
            posx_target = posx(posidx);
            posy_target = posy(posidx);
        end
        posidx = posidx+1;
        if posidx>8
            posidx=1;
        end
    end
    screen('flip',windowPtr,nextFlipTime);
    nextFlipTime = getSecs + settings.stimtime/1000;
    
    % SCREEN 3: DELAY
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
    screen('flip',windowPtr,nextFlipTime);
    pause(settings.delaytime/1000);
    
    % SCREEN 4: RESPONSE
    % Show a circle at the target location
    [mousex_bf mousey_bf buttons_bf] = GetMouse(windowPtr); 
    done_circle = 0;
    while ~done_circle
        screen('fillRect',windowPtr,settings.bgdac);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
        rect_target = CenterRectOnPoint([0 0 size(im)*1.2], posx_target, posy_target);
        screen('frameOval',windowPtr, 255, rect_target);
        screen('flip',windowPtr);

        [mousex_aft mousey_aft buttons_aft] = GetMouse(windowPtr);
        if (mousex_aft-mousex_bf)^2+(mousey_aft-mousey_bf)^2 > 25
            done_circle = 1;
        end
    end
    
    % response appears with mouse input
    done=0;
    mousex = screen_resolution(1)/2 + rand*screen_resolution(1)/4-screen_resolution(1)/8;        
    data.respstartangle(trialnr) = mod(mousex,screen_resolution(1)/4)/(screen_resolution(1)/4)*180-90;
    while ~done
        SetMouse(round(mousex),round(screen_resolution(2)/2),windowPtr);
        screen('fillRect',windowPtr,settings.bgdac);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
        % draw response gabor
        respangle = mod(mousex,screen_resolution(1)/4)/(screen_resolution(1)/4)*180-90;
        im = generate_gabor(respangle,settings.gablambda*screen_ppd,settings.gabsigma*screen_ppd,settings.gabphase);
        im = ((settings.gabpeak*im*256)+256)/2;
        patchsize = size(im);
        stimtex=screen('MakeTexture', windowPtr, im);
        srcrect = [0 0 patchsize];        
        destrect = centerRectOnPoint(srcrect,posx_target,posy_target);
        screen('drawtexture',windowPtr,stimtex,srcrect,destrect);               
%         screen('DrawText',windowPtr,['Angle= ' num2str(respangle,2)],0,0,[255 255 255]);
        screen('flip',windowPtr);
        
        [mousex mousey buttons] = GetMouse(windowPtr);
        done = any(buttons);
        if sum(buttons)>1
            screen('closeall');
            error('Experiment aborted');
        end
    end
    data.respangle(trialnr) = respangle;    

    % SCREEN 4: INTER TRIAL DISPLAY
    screen('fillRect',windowPtr,settings.bgdac);
    drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,0);
    screen('flip',windowPtr);
    nextFlipTime = getSecs + settings.intertrialtime/1000;
    
    if ~isempty(find(trialnr==breakpoints))
        save(outputfile,'settings','data');
        % BREAK
        breakStart=getSecs;
        while (GetSecs-breakStart)<settings.breakTime
            screen('fillRect',windowPtr,settings.bgdac);
            screen('DrawText',windowPtr,['You have finished ' num2str(round(100*trialnr/nTrials)) '% of this session.'],30,100,[255 255 255]);
            screen('DrawText',windowPtr,['Please take a short break now.'],30,120,[255 255 255]);
            totalBreak = getSecs-breakStart;
            screen('DrawText',windowPtr,['You can continue in ' num2str(ceil(settings.breakTime-totalBreak)) ' seconds.'],30,160,[255 255 255]);
            drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
            screen('flip',windowPtr);
        end
        screen('fillRect',windowPtr,settings.bgdac);
        screen('DrawText',windowPtr,['You have finished ' num2str(round(100*trialnr/nTrials)) '% of this session.'],30,100,[255 255 255]);
        screen('DrawText',windowPtr,['Please take a short break now.'],30,120,[255 255 255]);
        totalBreak = getSecs-breakStart;
        screen('DrawText',windowPtr,['Press any key to continue.'],30,160,[255 255 255]);
        screen('flip',windowPtr);
        waitForKey;
        screen('flip',windowPtr);
        drawfixation(windowPtr,screen_fixposxy(1),screen_fixposxy(2),250,5,1);
        screen('flip',windowPtr);        
        pause(1.5);
    end

end
save(outputfile,'settings','data');

%-%-%-%-%-%-%-
%- FINALIZE %-
%-%-%-%-%-%-%-
ShowCursor;
screen('closeall');



%-%-%-%-%-%-%-%-%-%-%-%-%- HELPER FUNCTIONS %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-

function keyCode = waitForKey
keyCode = ones(1,256);
while sum(keyCode(1:254))>0
    [keyIsDown,secs,keyCode] = KbCheck;
end
while sum(keyCode(1:254))==0
    [keyIsDown,secs,keyCode] = KbCheck;
end
keyCode = min(find(keyCode==1));

function drawfixation(windowPtr,x,y,color,size,vert)
Screen('DrawLine',windowPtr,color,x-size,y,x+size,y,2);
if (vert)
    Screen('DrawLine',windowPtr,color,x,y-size,x,y+size,2);
end

