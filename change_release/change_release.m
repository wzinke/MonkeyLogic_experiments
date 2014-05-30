% change_release (timing script)
%
% Initial timing file for teaching the NHP to use a joystick in order to get a reward.
%
% This timing file is prepared to run with the GitHub version of
% MonkeyLogic (https://github.com/Dfreedmanlab/MonkeyLogic_stable). It is possible that it will cause
% errors when used with previous versions of MonkeyLogic.
%
% This task utilizes two Stimuli. One of them signals the animal not to pull the joystick (NoGo Stimulus),
% whereas the other stimulus (Go Stimulus) requires a pull.
% The combination of both stimuli results in 3 task conditions:
%     1: Go only   - Only one stimulus is shown and the joystick needs to be pulled to get a reward
%     2: NoGo only - Only the NoGO stimulus is shown the animal has to withhold any responses to get a reward.
%     3: NoGo/Go   - First, the NoGO stimulus is shown, and after a variable time interval changes to the Go stimulus.
%                    The animal has to pull the joystick as response to the Go stimulus to receive a reward.
%
% The NoGo stimulus is shown for a random interval in the range of <dimmMin> and <dimmMax>.  The animal has to
% pull the joystick within a time window between <wait_resp> and <max_resp> following the onset of the Go stimulus
% to receive a reward.
%
% Be aware, that this timing file is made to be used with a analogue joystick as response device.
% The use of an analogue joystick changes the calls and logic of the eyejoytrack function.
%
% wolf zinke, May. 2014

%% ### ToDo: ### %%
    % use the gen function in the condition file to create the target item and make its features (size, luminance) editable.
    % define ITI as variable in TrialRecord (done!) and use the ITI end time instead of ML iti mechanism.
    % get subject name and paradigm from MLConfig (not accessible from timing file right now).

eventmarker(15);  %  start pre-trial

%% Define Trial Variables
tblpath = pwd;  % path to directory for trial table files
cNHP   = 'NHP_name';      % MLConfig.SubjectName    - MLConfig could not be accessed from the timing file
cPRDGM = 'get_joystick';  % MLConfig.ExperimentName - MLConfig could not be accessed from the timing file

% set editable variables
editable('Jpos0x', 'Jpos0y', 'JoyPullRad', 'JoyRelRad', 'dimmMin', 'dimmMax', 'dimmStep', 'wait_resp', 'max_resp', 'wait_rel', 'time_out', 'pull_pause', 'minITI', 'maxITI', 'ITIstep', 'RewInc2P', 'RewInc3P', 'RewInc4P', 'RewInc5P');

% initialize the random number generator
rng('shuffle', 'twister');

% reward times
rew_dur    =    50;    % duration of reward pulse (not to be changed!)
rew_gap    =   200;    % gap between two reward pulses (period = rew_dur+rew_gap) (not to be changed!)
rew_lag    =   100;    % introduce a brief delay prior reward
RewIncConsecutive = 1; % increase reward for subsequent correct trials. Otherwise reward will increase with the number of hits

% set centre position of joystick
Jpos0x     =   0;      % zero position X
Jpos0y     =   0;      % zero position Y
JoyPullRad = 1.5;      % threshold to detect a elevation of the joystick
JoyRelRad  =   1;      % threshold to detect the release of the joystick

% dimming parameters
dimmMin    =  250;     % minimum time to dimming
dimmMax    = 1500;     % maximum time to dimming
dimmStep   =   10;     % steps of possible dimming times

% inter trial times (replace this with a more sophisticated function)
minITI     =   500;    % minimum time period between two subsequent stimulus presentation
maxITI     =  3000;    % maximum time period between two subsequent stimulus presentation
ITIstep    =   100;

ITIvec  = minITI : ITIstep : maxITI;     % possible ITI's
cITI    = ITIvec(randi(length(ITIvec))); % ITI used after the current trial (uniform distribution)

% define ITI according to a shifted and truncated exponential distribution
%  mu = (minITI + maxITI)/2
%  max_range = exp(-(minITI/mu));
%  min_range = exp(-(maxITI/mu));
%  cITI = ITIstep * round( (-mu*reallog((max_range-min_range)*rand(1,1)+min_range)) / ITIstep);

% joystick response parameters
wait_resp  =  100;     % valid response only accepted after this initial period
max_resp   = 1500;     % maximal time accepted for a valid responses

time_out   =    0;     % set wait period for bad monkeys (do not combine with jittered ITIs)
pull_pause = 6000;     % this is the minimum time passed before a trial starts after random lever presses
wait_rel   = 1000;     % wait for release of joystick after response

%% Assign stimulus items
% assign items specified in the condition file to meaningful names
% numbers correspond to the number following TaskObject# in the condition file
JoyZero  =  1;         % this might be a simple trick to lock the joystick to a zero position by using a Fixation object
StimNoGo =  3;         % stimulus that indicates no response
StimGo   =  2;         % stimulus that requires a response

%% initialize trial variables
on_track   =   0;  % use this flag to indicate occurrences of errors in the trial
rt         = NaN;
RelTime    = NaN;
RewTime    = NaN;
WaitStart  = NaN;
TrialZero  = NaN;
TstartEff  = NaN;
TRelEff    = NaN;
DimmOnTime = NaN;
DimmEff    = NaN;
FixOn      = NaN;
FixOff     = NaN;

TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = NaN;
TrialRecord.Tend(  TrialRecord.CurrentTrialNumber) = NaN;

%% determine number of reward pulses (incremental reward scheme)
if(~isfield(TrialRecord, 'CorrCount'))
    TrialRecord.CorrCount = 0;
end

if(RewIncConsecutive == 1)
   % increase reward for consecutive correct trials
    RewInc2P = 1;   %   two pulses
    RewInc3P = 2;   % three pulses
    RewInc4P = 3;   %  four pulses
    RewInc5P = 4;   %  five pulses (jackpot)

    if(TrialRecord.CorrCount < RewInc2P-1)
        rew_Npulse = 1;
    elseif(TrialRecord.CorrCount < RewInc3P-1)
        rew_Npulse = 2;
    elseif(TrialRecord.CorrCount < RewInc4P-1)
        rew_Npulse = 3;
    elseif(TrialRecord.CorrCount < RewInc5P-1)
        rew_Npulse = 4;
    else
        rew_Npulse = 5;
    end
else
    cNumHit = sum(TrialRecord.TrialErrors == 0);

    % increase rewards after a defined number of trials was achieved
    RewInc2P =  75;   %   two pulses
    RewInc3P = 200;   % three pulses
    RewInc4P = 350;   %  four pulses
    RewInc5P = 450;   %  five pulses (jackpot)

    if(cNumHit < RewInc2P)
        rew_Npulse = 1;
    elseif(cNumHit < RewInc3P)
        rew_Npulse = 2;
    elseif(cNumHit < RewInc4P)
        rew_Npulse = 3;
    elseif(cNumHit < RewInc5P)
        rew_Npulse = 4;
    else
        rew_Npulse = 5;
    end
end

% split the stimulus change/dimming times across conditions to allow for some control
num_cnd = length(unique(TrialRecord.ConditionsThisBlock)); % This variable is behaving strange and not as I understand its role. The use of unique ensures that it shows the actual number of conditions used in this block!
ccnd    = 1+mod(TrialRecord.CurrentCondition-1,3) == 0;    % task conditions 1:3
dcnd    = ceil(TrialRecord.CurrentCondition/3);            % change duration block

if(ccnd > 1) % The first conditions just shows the Go-Stimulus
    itvtm   = linspace(dimmMin, dimmMax, (num_cnd/3)+1);
    Dimmvec = dimmMin : dimmStep : dimmMax; % possible stimulus change times
    Dimmvec = Dimmvec(Dimmvec >= itvtm(dcnd) & Dimmvec < itvtm(dcnd+1)); % possible stimulus change times for current condition

    cDimm   = Dimmvec(randi(length(Dimmvec)));  % stimulus change time used for the current trial
else
    cDimm = 0;
end

%  mu = (dimmMin + dimmMax)/2
%  max_range = exp(-(dimmMin/mu));
%  min_range = exp(-(dimmMax/mu));
%  R = dimmStep * round( (-mu*reallog((max_range-min_range)*rand(1,1)+min_range)) / dimmStep);

%%%%###########################################%%%%
%%%%########  Ensure joystick release   #######%%%%
%%%%###########################################%%%%
%% check for joystick position.
% If the joystick is pulled, a defined period will be waited for a release,
% and this trial is aborted either after this period passed or the lever is released.
% This trial restart is necessary to avoid ML crashes due too long trial durations.

reposition_object(JoyZero, Jpos0x, Jpos0y);  % correct for offsets of the joystick analogue output

if(~isfield(TrialRecord, 'Jreleased'))
    TrialRecord.Jreleased = 1;   % track the joystick state across trials
end

% check if lever remains released
[JoyRel PullTime] = eyejoytrack('holdtarget', JoyZero, JoyRelRad, 100); % this introduces a 100 ms interval

if(~JoyRel | TrialRecord.Jreleased == 0)
    if(~JoyRel) % wait for release
        TrialRecord.Jreleased = 0;
        [JoyRel] = eyejoytrack('acquiretarget', JoyZero, JoyRelRad, pull_pause);

        if(JoyRel)
            eventmarker(39); % lever released
        end

    else % make sure that joystick is not pulled for a sufficient long time before starting a new trial
        [JoyRel PullTime] = eyejoytrack('holdtarget', JoyZero, JoyRelRad, pull_pause);

        if(JoyRel)
            TrialRecord.Jreleased = 1;
        end
    end

    TrialRecord.CorrCount = 0;
    trialerror(8);
    set_iti(0);
else
    on_track = 1;
    TrialRecord.Jreleased = 1;
end

eventmarker(16); % End pre trial

%%%%###########################################%%%%
%%%%##########  Start the Trial now   #########%%%%
%%%%###########################################%%%%

if(on_track)
%% #### start Trials #### %%
    eventmarker(10);  % end ITI

    if(~isfield(TrialRecord, 'NextTrial'))
        TrialRecord.NextTrial = -1;
    end

    if(ccnd > 1) % The first conditions just shows the Go-Stimulus, so skip this condition here
%% #### Onset NoGo Stimulus #### %%
        NoGoOn = toggleobject(StimNoGo, 'EventMarker', [1, 35], 'Status', 'on');

        % get the computer time for the trial start as well
        TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = toc(uint64(1));
        TstartEff = TrialRecord.Tstart(TrialRecord.CurrentTrialNumber);

        FixOn     = NoGoOn;
        TrialZero = NoGoOn;

%% #### Wait Change #### %%
        [JoyHold ReleaseTime] = eyejoytrack('holdtarget', JoyZero, JoyPullRad, cDimm);  % wait before stim change (hope no time passed since onset)

        % do not accept early responses
        if(~JoyHold)
            eventmarker(38);  % lever pressed

            % get the computer time for the lever release as well
            TPullEff = toc(uint64(1));
            RelTime =  trialtime - TrialZero;

            TrialRecord.Jreleased = 0;
            ErrCode  = 5;
            on_track = 0;
            StimOffTime = toggleobject(StimNoGo , 'EventMarker', 36,'Status', 'off');
        end
    end  % if(ccnd > 1)

%% #### Onset Go Stimulus #### %%
    % wait for release of joystick
    if(on_track & ccnd ~= 2) % skip the NoGo trial
        if(ccnd == 1)
            GoOn = toggleobject([StimNoGo StimGo], 'EventMarker', [37]);  % show Go stimulus
            StartChngWait = toc(uint64(1));
        else
            GoOn = toggleobject(StimGo , 'EventMarker', [37, 35, 1],'Status', 'on');
            StartChngWait = toc(uint64(1));

            % get the computer time for the trial start as well
            TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) = toc(uint64(1));
            TstartEff = TrialRecord.Tstart(TrialRecord.CurrentTrialNumber);
            FixOn     = GoOn;
            TrialZero = FixOn;
        end

        [JoyHold PullTime] = eyejoytrack('holdtarget', JoyZero, JoyRelRad, max_resp);
        EndChngWait = toc(uint64(1));

        if(JoyHold)
            ErrCode  = 1; % no response
            on_track = 0;
            StimOffTime = toggleobject(StimGo, 'EventMarker', 36,'Status', 'off');
        else
            eventmarker(38);  % lever pressed

            PullTime = trialtime - TrialZero;

            TPullEff = EndChngWait;

            TrialRecord.Jreleased = 0;

            if(PullTime < wait_resp)
                ErrCode  = 5; % early response
                on_track = 0;
                StimOffTime = toggleobject(StimGo , 'EventMarker', 36,'Status', 'off');
            else
                rt = PullTime % double check this one (use text table and compare to RelTime-DimmOnTime!)
                rt = 1000 * (EndChngWait - StartChngWait) % double check this one (use text table and compare to RelTime-DimmOnTime!)
            end
        end
    end  % if(on_track & ccnd ~= 2)


%% #### Trial End #### %%
   eventmarker([2 17]);  % end of trial

    % finalize Trial
    if(on_track)  % correct trials
        ErrCode = 0;
        trialerror(0);  % correct

        TrialRecord.CorrCount = TrialRecord.CorrCount + 1;

        eventmarker(19);   % start pause
        idle(rew_lag);

        if(ccnd == 2)
            StimOffTime = toggleobject(StimNoGo, 'EventMarker', 36,'Status', 'off');
        else
            StimOffTime = toggleobject(StimGo, 'EventMarker', 36,'Status', 'off');
        end
        eventmarker([96 20]);

        RewEff = toc(uint64(1));
        goodmonkey(rew_dur, 'NumReward', rew_Npulse, 'PauseTime', rew_gap);
        eventmarker(96); % use this event twice to get an estimate of the reward duration, i.e. number of pulses.
        RewTime = trialtime - TrialZero;

        if(TrialRecord.Jreleased == 0)
            [JoyRel] = eyejoytrack('acquiretarget', JoyZero, JoyRelRad, wait_rel);
            if(JoyRel)
                eventmarker(39); % lever released
                EndChngWait = toc(uint64(1));
                TrialRecord.Jreleased == 1;
            end
        end
    else   % error occurred
        trialerror(ErrCode);
        rew_Npulse = NaN;
        TrialRecord.CorrCount = 0;
        cITI = cITI + time_out;   % add a punish time to the current ITI
    end  %  if(on_track)

    FixOff = StimOffTime - TrialZero;

    % get the computer time for the trial end as well
    TrialRecord.Tend(TrialRecord.CurrentTrialNumber) = toc(uint64(1));

%%%%###########################################%%%%
%%%%##########   Stop the Trial now   #########%%%%
%%%%###########################################%%%%
    % set ITI
    eventmarker([18 9]); % end post trial / start ITI

    if(TrialRecord.CurrentTrialNumber > 1)
        ITIeff = TrialRecord.Tstart(TrialRecord.CurrentTrialNumber) - TrialRecord.Tend(TrialRecord.CurrentTrialNumber-1);
    else
        ITIeff = NaN;
    end

    set_iti(cITI);
    TrialRecord.NextTrial = TrialRecord.Tend(TrialRecord.CurrentTrialNumber) + cITI/1000;

    if(ErrCode ~= 8 && ErrCode ~= 9)
    %% write trial table
    % keep redundant information to check data and have some information as back up to reconstruct trials
    % Further, this text file should give fast and easy access to the behavioral performance. It is in a
    % format that allows to be read into R with the read.table command (use 'header=TRUE' option).
        cdt   = datestr(now,'yyyy_mm_dd');
        tblnm = fullfile(tblpath, [cNHP,'_', cPRDGM,'_', cdt, '.dat']);

        if(exist(tblnm) ~= 2 )  % create the table file
            tblptr = fopen(tblnm, 'w');
            fprintf(tblptr,'Date  Subject  Experiment  TFixOn  Tstart  TRel  Tend  ITIeff  TrialNo  BlockNo  CondNo  Result  WaitStart  ITI  NumRew  RT  TrialZero  FixOn  FixOff  DimmEff  DimmOn  DimmIntend  ReleaseTime  RewTime\n');
        else
            tblptr = fopen(tblnm, 'a');
        end

        fprintf(tblptr,'%s  %s  %s  %.4f  %.4f  %.4f  %.4f  %6d  %6d  %4d  %4d  %4d  %8d  %8d  %8d  %8d  %8d  %8d  %8d  %.4f  %8d  %8d  %6d  %6d \n',...
                cdt, cNHP, cPRDGM, TrialRecord.Tstart(TrialRecord.CurrentTrialNumber), TstartEff, TRelEff, ...
                TrialRecord.Tend(TrialRecord.CurrentTrialNumber), ITIeff, TrialRecord.CurrentTrialNumber,  ...
                TrialRecord.CurrentBlock, ccnd, ErrCode, WaitStart, cITI, rew_Npulse, ...
                rt, TrialZero, FixOn, FixOff, DimmEff, DimmOnTime, cDimm, RelTime, RewTime );
        fclose(tblptr);
    end
end

