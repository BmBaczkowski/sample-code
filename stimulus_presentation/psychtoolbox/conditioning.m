function conditioning(varargin)

%% STEP 0: CONFIG PARAMETERS

% Clear the screen
sca;
clc;

% Print info to the command window
more off;

% set default parameters
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
screenNumber = 1;
debugging = 1;
debugging_screen = 0;
dummymode_tracker = 1;
dummymode_pp = 1;
inputkeyboard = 'Dell Dell Wired Multimedia Keyboard';
%inputkeyboard = 'HID 046a:0023';
%inputkeyboard = 'Teensyduino MilliKey';
%inputkeyboard = 'Magic Keyboard';
%datdir changed for github situation
datdir = 'data/temp';
basedir = './';

pp = nan;
session_name = 'conditioning';

% read in the parameters from the function arguments
for i=1:length(varargin)
    switch cell2mat(varargin(i))
        case 'dummymode_tracker'
            dummymode_tracker = varargin{i+1};
        case 'dummymode_pp'
            dummymode_pp = varargin{i+1};
        case 'screenNumber'
            screenNumber = varargin{i+1};
        case 'datdir'
            datdir = varargin{i+1};
        case 'inputkeyboard'
            inputkeyboard = varargin{i+1};
        case 'debugging_screen'
            debugging_screen = varargin{i+1};
        case 'debugging'
            debugging = varargin{i+1};
        case 'session_name'
            session_name = varargin{i+1};
        case 'basedir'
            basedir = varargin{i+1};
    end
end

resourcesdir = fullfile(basedir, 'resources');
stimsdir  = fullfile(basedir, 'stims');

% load helper functions
addpath(genpath(resourcesdir));

% setup the input keyboard for reading in responses
KbPointer = GetKeyboardIndices(inputkeyboard);
KbPointer = KbPointer(1);

%% STEP 1: LOGFILES

% get participant id
prompt = {'Enter participant number (3 characters)'};
dlg_title = 'Welcome';
def = {'000'};
answer = inputdlg(prompt, dlg_title, 1, def);

% Print some text in Matlab's Command Window if a file name has not been entered
if  isempty(answer)
    fprintf('Session cancelled by user\n')
    cleanup(dummymode_pp, pp);
    return
end

% read in subject id
subjectid = answer{1};
subjectdir = fullfile(datdir, subjectid);

% get subject directory
if ~(exist(subjectdir, 'dir') == 7)
    mkdir(subjectdir);
end

% logfile
logfilename     =   sprintf('%s.csv', session_name);
logfile         =   fopen(fullfile(subjectdir, logfilename),'w+');
fprintf(logfile, '%s,', ...
    'TRIAL', ...
    'PICTURE', ...
    'NODE', ...
    'GRAPH', ...
    'CS', ... % 1-CS+; 2-CS-
    'SHOCK', ... % 0-no shock; 1-shock
    'MARKER', ...
    'RESPONSE', ... %  (1-expected shock; 0-did not)
    'RT');
fprintf(logfile, '\n');
fclose(logfile);

%% STEP 2: RANDOMIZATION

load(fullfile(resourcesdir, 'pseudorandTrials4conditioning.mat'), 'pseudorandTrials');
% select randomly a set of trials
trials_conditions   = pseudorandTrials(:,randi(size(pseudorandTrials,2)));
% make the two column var
conditions              = [trials_conditions zeros(size(trials_conditions,1),1)]; %column 1 CS+/CS-; column2 shock/noshock
conditions(conditions(:,1)==3,1)          = 1;
conditions(trials_conditions==3,2)    = 1;
clear pseudorandTrials trials_conditions;

if debugging
    nTrials = 2; % number has to be even
else
    nTrials = size(conditions,1);
end

% Get 'design mat' n x 6 [trial, picture, node, cs, shock, marker]
X               = nan(nTrials,6);

% prepare inter trial interval
if isOctave
    pkg load statistics;
end

ITI = nan(nTrials,1);
while nanmean(ITI) ~= 9 % mean of all ITIs is 9 sec
    for i=1:nTrials
        ITI(i) = randsample([7,9,11],1);
    end
end

% select CS+ and CS- nodes
nodes = [1, 13];
% select graph
graph = 1;
% read pics2node allocation
pics2nodes_filename = fullfile(subjectdir, 'pics2nodesObj.json');
if exist(pics2nodes_filename, 'file') == 2
    pics2nodes = get_pics2nodes(pics2nodes_filename);
else
    fprintf('Pics2nodes file does not exist.\nSession cancelled.\n');
    cleanup(dummymode_pp, pp);
    return
end

%% STEP 3: INIT PARALLEL PORT & EYELINK

% Init parallel port
[dummymode_pp, pp] = initParallelPort('dummymode', dummymode_pp);
if dummymode_pp
    fprintf('Running parallel port in a dummy mode\n\n\n');
end

% Init Eyelink
EyelinkInit(dummymode_tracker);

% get filename for eyelink edf
edfFile = sprintf('%s', session_name(1:3));

% Open an EDF file and name it
failOpen = Eyelink('OpenFile', edfFile);
if failOpen ~= 0 % Abort if it fails to open
    fprintf('Cannot create EDF file ''%s', edfFile);
    cleanup(dummymode_pp, pp);
    return
end

% Get EyeLink tracker and software version
% <ver> returns 0 if not connected, returns 1 for EyeLink I, 2 for EyeLink II, 3 for EyeLink 1K/Plus
% <versionstring> returns 'EYELINK I', 'EYELINK II x.xx', 'EYELINK CL x.xx' where 'x.xx' is the software version
[ver, versionstring] = Eyelink('GetTrackerVersion');
if ~dummymode_tracker
    fprintf('Running experiment on %s version %d\n', versionstring, ver );
end

% Add a line of text in the EDF file to identify the session. This is optional
preambleText = sprintf('PSYCHTOOLBOX EXPERIMENT: %s', session_name);
Eyelink('Command', 'add_file_preamble_text "%s"', preambleText);
% Select which events are saved in the EDF file. Include everything just in case
Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
% Select which events are available online for gaze-contingent experiments. Include everything just in case
Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,FIXUPDATE,INPUT');
% Select which sample data is saved in EDF file or available online. Include everything just in case
if ver > 2  % Check tracker version and include 'HTARGET' to save head target sticker data for supported eye trackers
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,HTARGET,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
else
    Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,RAW,AREA,GAZERES,BUTTON,STATUS,INPUT');
    Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT');
end
% set the sample rate
Eyelink('Command', 'sample_rate = 500');

%% STEP 4: OPEN PTB WINDOW

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);

% background color
grey = 128/255;

% Open the window
if debugging_screen
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [0 0 800 600]);
else
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey);
end

% clear the screen
Screen('Flip', window);

% makes it so that characters typed don't show up in the command window
ListenChar(-1);

% hide cursor
HideCursor(window);

% Set up alpha-blending for smooth (anti-aliased) lines to draw fixation cross
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% GET WINDOW PROPERTIES

% Return width and height of the graphics window/screen in pixels
[width, height] = Screen('WindowSize', window);
%Retrieve monitor refresh rate
hz = Screen('NominalFrameRate', window);
% Get the centre coordinate of the window in pixels.
[xCenter, yCenter] = RectCenter(windowRect);

% Query the inter-frame-interval. This refers to the minimum possible time
% between drawing to the screen
ifi = Screen('GetFlipInterval', window);

% Get the maximum coded luminance level (this should be 1)
maxLum = Screen('ColorRange', window);

% Print resolution and refreh rate in Matlab's Command Window
WaitSecs(1);
fprintf('\n\n\n\n');
fprintf('Running on a %d x %d screen at %d Hz\n', width, height, hz);
fprintf('Inter-frame-interval is %d\n', ifi);
fprintf('Maximum coded luminance is %d\n', maxLum);
fprintf('\n\n\n\n');

%% STEP 5: LOAD IMAGES INTO TEXTURES

% Set the text format
Screen('TextSize', window, 42);
Screen('TextFont', window, 'Arial');
Screen('TextStyle', window, 1); % Bold

dotColor    = [68/255 68/255 248/255];
dotXpos     = xCenter;
dotYpos     = yCenter;
%Fixation dot of .7 degrees visual angle is 26 pixels.
dotSizePix  = 26;

% load images
stimTex = nan(size(nodes));
for i = 1:size(nodes,2)
    pic           = pics2nodes(nodes(i),graph);
    imagefilename = fullfile(stimsdir, sprintf('hsv_n_%02d.png', pic));
    img           = imread(imagefilename);
    imgSize       = size(img,1);
    stimTex(i)    = Screen('MakeTexture', window, img);
end
clear img;

%% STEP 6: SHOW WELCOME SCREEN AND INSTRUCTIONS

text = 'WELCOME';
DrawFormattedText(window, text, 'center', 'center', .8);
Screen('Flip', window);
KbStrokeWait();


%% STEP 7: SET EYELINK CALIBRATION

% Provide EyeLink with some defaults, which are returned in the structure "el".
el = EyelinkInitDefaults(window);
% set calibration/validation background and target colors
el.backgroundcolour = grey; % RGB grey
el.calibrationtargetcolour = 0; % RGB black
% Set calibration beeps (0 = sound off, 1 = sound on)
el.targetbeep = 1;    % sound a beep when a target is presented
el.feedbackbeep = 1;  % sound a beep after calibration or drift check/correction
% You must call this function to apply the changes made to the el structure above
EyelinkUpdateDefaults(el);

% Set display coordinates for EyeLink data by entering left, top, right and bottom coordinates in screen pixels
Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
% Write DISPLAY_COORDS message to EDF file: sets display coordinates in DataViewer
% See DataViewer manual section: Protocol for EyeLink Data to Viewer Integration > Pre-trial Message Commands
Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);

% Set number of calibration/validation dots and spread: horizontal-only(H) or horizontal-vertical(HV) as H3, HV3, HV5, HV9 or HV13
Eyelink('Command', 'calibration_type = HV9'); % horizontal-vertical 9-points
% Clear Host PC display from any previus drawing
Eyelink('Command', 'clear_screen 0');
% Put EyeLink Host PC in Camera Setup mode for participant setup/calibration
EyelinkDoTrackerSetup(el);


%% STEP 8: SHOW INFO TO THE SUBJECT (and start the task)

% Show instructions
text = 'The task will start shortly...';
DrawFormattedText(window, text, 'center', 'center', .8);
Screen('Flip', window);
fprintf('Press SPACE to start the task\n\n');
KbStrokeWait();
Screen('Flip', window);


%% STEP 9: EYELINK PREP

% Put tracker in idle/offline mode before recording
Eyelink('SetOfflineMode');
% Clear Host PC display from any previus drawing
Eyelink('Command', 'clear_screen 0');
% draw a box on the host PC to indicate the image outline
Eyelink('Command', 'draw_box %d %d %d %d 15', round(width/2-imgSize/2), round(height/2-imgSize/2), round(width/2+imgSize/2), round(height/2+imgSize/2));
% draw a cross in the center of the image
Eyelink('Command', 'draw_cross %d %d 15', width/2 - 1, height/2 - 1);
% Allow some time for drawing
WaitSecs(0.1);

% Start tracker recording
Eyelink('StartRecording');
% Allow some time
WaitSecs(0.1);
% Check the recording status
status = Eyelink('CheckRecording');
if status~=0
    fprintf('Eyelink is not recording.\n');
    cleanup(dummymode_pp, pp);
    return;
end

%% STEP 10: START THE TASK
Screen('Flip', window);

% record a few samples before we actually start displaying
% mark zero-plot time in data file
Eyelink('Message', 'SYNCTIME');
WaitSecs(1);

% Mark the begining
KbQueueCreate(KbPointer);
KbQueueStart(KbPointer);
Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
Screen('Flip', window);
sendMarker(dummymode_pp, pp, 101);
Eyelink('Message', '101');
t = clock;
fprintf('TASK START at %02d:%02d\n\n', t(4), t(5));
WaitSecs(10);

%% STEP 11: TRIAL LOOP

for trial = 1:nTrials
    KbEventFlush(KbPointer);
    
    counter         = trial;
    X(counter,1)    = trial;                                          % which trial
    X(counter,2)    = pics2nodes(nodes(conditions(trial,1)), graph);  % which picture
    X(counter,3)    = nodes(conditions(trial,1));                     % which node
    X(counter,4)    = graph;                                          % which graph
    X(counter,5)    = conditions(trial,1);                            % which CS
    X(counter,6)    = conditions(trial,2);                            % shock
    X(counter,7)    = (X(counter,5) * 10) + (X(counter,6) * 20); % marker 10- CS+; 20- CS-; 30- CS+UCS+
    
    % Draw stimulus
    Screen('DrawTexture', window, stimTex( X(counter,5) )); % texture is the same as CS
    % Draw fixation dot
    Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
    KbQueueFlush(KbPointer);
    stimOnset = Screen('Flip', window);
    sendMarker(dummymode_pp, pp, X(counter,7));
    Eyelink('Message', num2str(X(counter,7)));
    
    Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
    Screen('DrawingFinished', window);
    
    if X(counter,6) == 1
        WaitSecs('UntilTime', stimOnset+3.8);
        Eyelink('Message', '128');
        sendShock(dummymode_pp, pp);
    end
    
    stimOffset = Screen('Flip', window, stimOnset+4.0-ifi/2);
    
    [ keyIsDown, keyCode ] = KbQueueCheck(KbPointer);
    % find the first button that was pressed
    keyCode(keyCode==0) = inf; 
    keyCodeIndx = find(keyCode==min(keyCode));
    response = NaN;
    RT = NaN;
    if keyIsDown
        RT = keyCode(keyCodeIndx) - stimOnset;
        if keyCodeIndx == KbName('1!')
            response = 1;
        elseif keyCodeIndx == KbName('3#')
            response = 0;
        else
            RT = NaN;
            fprintf('Trial %d: Wrong button was pressed!!!\n\n', trial);
        end
    else
        fprintf('Trial %d: No button press!!!\n\n', trial);
    end
    
    logfile =  fopen(fullfile(subjectdir, logfilename),'a+');
    fprintf(logfile, '%s,', ...
        num2str(X(counter,1)), ...
        num2str(X(counter,2)), ...
        num2str(X(counter,3)), ...
        num2str(X(counter,4)), ...
        num2str(X(counter,5)), ...
        num2str(X(counter,6)), ...
        num2str(X(counter,7)), ...
        num2str(response), ...
        sprintf('%0.3f', RT));
    fprintf(logfile, '\n');
    fclose(logfile);
    
    WaitSecs('UntilTime', stimOffset+ITI(trial));
    
    % add break in between
    if (debugging == 1 && trial == nTrials/2) || (debugging == 0 && trial == nTrials/2)
        
        text = 'Pause\n\n (12 Sekunden)';
        DrawFormattedText(window, text, 'center', 'center', .8);
        Screen('Flip', window);
        fprintf('Break (12 seconds)\n');
        sendMarker(dummymode_pp, pp, 110);
        Eyelink('Message', '110');
        Screen('DrawDots', window, [dotXpos dotYpos], dotSizePix, dotColor, [], 2);
        Screen('DrawingFinished', window);
        WaitSecs(12);
        fprintf('Break end\n\n');
        Screen('Flip', window);
        sendMarker(dummymode_pp, pp, 111);
        Eyelink('Message', '111');
        WaitSecs(6);
    end
end
KbQueueStop(KbPointer);
WaitSecs(10);
sendMarker(dummymode_pp, pp, 102);
Eyelink('Message', '102');



%% STEP 12: STOP EYELINK RECORDING, CLOSE & CLEAN
% Stop tracker recording
Eyelink('StopRecording');
% Allow some time
WaitSecs(0.1);

% Put tracker in idle/offline mode
Eyelink('SetOfflineMode');
% Close EDF file on Host PC
Eyelink('CloseFile');
% Clear trial image on Host PC at the end of the experiment
Eyelink('Command', 'clear_screen 0');
% Allow some time for screen drawing
WaitSecs(0.1);

% Transfer a copy of the EDF file to Display PC
transferFile(edfFile, subjectdir);
WaitSecs(.1);
cleanup(dummymode_pp, pp);

end