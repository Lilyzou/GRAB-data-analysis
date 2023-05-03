% This file will run automated; please copy stdshade.m into the directory
% first; uncomment standard error section and resave as semshade.m
% use 0969_0115_001 as an example

%10/15 notes: This version includes base, 16 frame normalization, and 200
%ms before stimuli normalization. Includes average responses of Hit, Miss,
%FA, and CR trials, with error. Includes movement analysis. Aligned to pole onset, pole down, first
%lick, first reward lick, last lick.

%10/25 notes: Fixes axes and adds true first lick (not reward)

%% load all the files that has info you need
tic

disp('Enter each session as the system prompts. To finish, type "go"');
warning('off','all');
stopRead = false;
i = 1;
fileNames(1) = strings;
while (stopRead == false)
    fNames = input(sprintf('Session number %d:', i), 's');
    fileNames(i) = fNames;
    if (strcmp(fileNames(i),'go'))
        stopRead = true;
        fileNames = fileNames(1:(i-1));
    else
        i = i+1;
    end
end

timeNormData = cell(length(fileNames),23); %23 represents the number of hm arrays that are generated in each loop
frameNormData = cell(length(fileNames),23);
baseData = cell(length(fileNames),23);

for sessionLoop = 1:length(fileNames)
    fname = char(fileNames(sessionLoop));
    disp(fname);
    fmat = strcat(fname,'.mat');
    ftrial = strcat(fname,'.trials');
    behave = fname(5:9);
    behave = strcat('behavior',behave,'.mat');
    fAlign = strcat(fname,'.align');

    currFolder = pwd;
    cd(strcat(currFolder,'\',fname(6:9),'\'));
    
    mkdir('trial-data-mean-base');
    mkdir('trial-data-mean-16frames');
    mkdir('trial-data-mean-200ms');
    mkdir('movement-analysis');
    mkdir('movement-analysis-16frames');
    mkdir('movement-analysis-200ms');
    currFolder = pwd;
    baseMeanFolder = strcat(currFolder,'\trial-data-mean-base\');
    normMeanFolder = strcat(currFolder,'\trial-data-mean-16frames\');
    timeMeanFolder = strcat(currFolder,'\trial-data-mean-200ms\');
    moveFolder1 = strcat(currFolder,'\movement-analysis\');
    moveFolder2 = strcat(currFolder,'\movement-analysis-16frames\');
    moveFolder3 = strcat(currFolder,'\movement-analysis-200ms\');
    whiskname = strcat('AH',fname(1:4),fname(6:9),'-WTLIA.mat');
    loadname = strcat(fname(6:9),'load.mat');
    if isfile(loadname) %%%%%%%%%%%%%%%IF YOU CHANGE ANYTHING WITHIN THIS IF ELSE BLOCK, MAKE SURE TO DELETE THE LOAD FILE SO THAT YOU CAN GENERATE A NEW FILE
        % if exist(loadname, 'file') == 2   %for matlab2015
        load (loadname);
        currFolder = pwd;
        baseMeanFolder = strcat(currFolder,'\trial-data-mean-base\');
        normMeanFolder = strcat(currFolder,'\trial-data-mean-16frames\');
        timeMeanFolder = strcat(currFolder,'\trial-data-mean-200ms\');
        moveFolder1 = strcat(currFolder,'\movement-analysis\');
        moveFolder2 = strcat(currFolder,'\movement-analysis-16frames\');
        moveFolder3 = strcat(currFolder,'\movement-analysis-200ms\');
    else
        %% get start frame and end frame of each trial from trial file
        %find the start_idx and end_idx
        disp('loading files')
        load(fmat)
        load(ftrial,'-mat')
        load(behave)
        if isfile(fAlign)
            load(fAlign, '-mat')
        end
        for i =1:length(trials)
            start_idxog(i)=trials(i).frames(1);
            end_idxog(i)=trials(i).frames(2);
        end
        %get the trial length for each trial
        for i=1:length(start_idxog)
            trial_lengthog(i)=end_idxog(i)-start_idxog(i);
        end
        %% pull out information from the raw file
        %can either read the raw file or the aligned file(_rigid.sbx)

        %automate file read
        %x=sbxread(strcat(fname,'_rigid'),0,end_frame);
        %y=sbxread(strcat(fname,'_rigid'),0,200);

        %best way:
        if isfile(strcat(fname,'_rigid.mat'))
            disp('loading rigid')
            rigidName = strcat(fname,'_rigid');
            readLoops = floor(end_frame/10000);
            startRead = 0;
            x1=sbxread(rigidName,startRead,10000);
            z1=squeeze(x1(1,:,:,:));%(choose the ROI at 2nd and 3rd dimension)
            tmp1=squeeze(mean(mean(z1,1),2));
            clear('x1','z1')
            startRead = startRead+10000;

            fprintf('load is %.2f%% complete\n',1/(readLoops+1)*100);
            for (i = 2:readLoops)
                x1=sbxread(rigidName,startRead,10000);
                z1=squeeze(x1(1,:,:,:));%(choose the ROI at 2nd and 3rd dimension)
                tmp2=squeeze(mean(mean(z1,1),2));
                clear('x1','z1')
                tmp1 = cat(1,tmp1,tmp2);
                startRead = startRead+10000;
                fprintf('load is %.2f%% complete\n',i/(readLoops+1)*100);
            end
            x1 = sbxread(rigidName,startRead,end_frame-startRead);
            z1 = squeeze(x1(1,:,:,:));
            tmp2 = squeeze(mean(mean(z1,1),2));
            clear('x1','z1')
            tmp = cat(1,tmp1,tmp2);
            disp('load is 100% complete')

            y=sbxread(strcat(fname,'_rigid'),0,200);
            save (loadname, '-regexp','^(?!(fileNames|timeNormData|baseData|frameNormData)$).');

        else
            disp('loading mat')
            readLoops = floor(end_frame/10000);
            startRead = 0;
            x1=sbxread(fname,startRead,10000);
            z1=squeeze(x1(1,:,:,:));%(choose the ROI at 2nd and 3rd dimension)
            tmp1=squeeze(mean(mean(z1,1),2));
            clear('x1','z1')
            startRead = startRead+10000;

            fprintf('load is %.2f%% complete\n',1/(readLoops+1)*100);
            for (i = 2:readLoops)
                x1=sbxread(fname,startRead,10000);
                z1=squeeze(x1(1,:,:,:));%(choose the ROI at 2nd and 3rd dimension)
                tmp2=squeeze(mean(mean(z1,1),2));
                clear('x1','z1')
                tmp1 = cat(1,tmp1,tmp2);
                startRead = startRead+10000;
                fprintf('load is %.2f%% complete\n',i/(readLoops+1)*100);
            end
            x1 = sbxread(fname,startRead,end_frame-startRead);
            z1 = squeeze(x1(1,:,:,:));
            tmp2 = squeeze(mean(mean(z1,1),2));
            clear('x1','z1')
            tmp = cat(1,tmp1,tmp2);
            disp('load is 100% complete')

            y=sbxread(fname,0,200);
            save (loadname, '-regexp','^(?!(fileNames|timeNormData|baseData|frameNormData)$).');
        end



    end
    toc
    %% Visualize data, photobleach elimination, FI trace
    n = 1; %%%%%%%%%%%%%%NOTE N HERE DENOTES FIGURE NUMBER

    %Automate photobleach elimination of trials; current threshold set to 100
    %FI gap, but this can be adjusted as needed
    % delTrial = 0;
    % for i = 1:30
    %     if (max(tmp(start_idxog(i):end_idxog(i)))-min(tmp(start_idxog(i+1):end_idxog(i+1))) > 50)
    %         delTrial = i;
    %     end
    % end
    % delTrial = delTrial + 1;
    % delTrial
    % text(end_idxog(delTrial),tmp(end_idxog(delTrial))*1.25,strcat(num2str(delTrial),' trials deleted'));
    delTrial = 30; %just for now


    start_idx=start_idxog(delTrial:end-1);
    end_idx=end_idxog(delTrial:end-1);
    trial_length=trial_lengthog(delTrial:end-1);
    %get F.I of each trial
    tmp2 = cell(length(start_idx),1);
    for i = 1:length(start_idx)
        tmp2{i} = tmp(start_idx(i):end_idx(i));
    end
    %plot the df trace of each trial

    %% plot the heatmap
    maxTmpLength = 0;
    for i = 1:length(start_idx)
        if length(tmp2{i}) > maxTmpLength
            maxTmpLength = length(tmp2{i});
        end
    end
    hm = NaN(length(start_idx),maxTmpLength);
    tmp3 = NaN(length(tmp2));
    for i = 1:length(start_idx)
        hm(i,1:length(tmp2{i})) = tmp2{i}/mean(tmp2{i})-1;
        tmp3(i) = mean(tmp2{i});
    end
    hm2=hm(:,1:100);
    frameRate = 15.44;
    secResolution = 1;


    % set the time from 1-5s
    frameRate = 15.44;
    secResolution = 1;

    %use the first 16 frames as baseline
    hm3 = NaN(length(start_idx),maxTmpLength);
    for i = 1:length(start_idx)
        hm3(i,1:length(tmp2{i})) = (tmp2{i}/mean(tmp2{i}(1:16)))-1;
    end

    hm4 = hm3(:,1:100);

    %pull out the first 100 frames from each normalized trial
    %generate the heatmap of first 100 frames
 


    %% plot the sorted heatmap by trial type
    %Automated Sorting by trial type with average responses ready to be plotted

    %Hit Trials
    hmbsort = zeros(length(hm2)+200,100); % generate an empty array to hold the sorted values
    counter1 = 1;
    counter2 = 1;
    hitSize = length(b{1}.hitTrialNums);
    hitCount = 0;
    if ~(isempty(b{1}.hitTrialNums))
        while (b{1}.hitTrialNums(end-hitCount)-delTrial >= length(start_idx)) % remove last trial if it trails off and doesn't end
            hitSize = hitSize-1;
            hitCount = hitCount+1;
        end
    end
    % initialize variables of start points and trials to delete
    hitDel = 1;
    missDel = 1;
    crDel = 1;
    faDel = 1;
    missStart = 1;
    crStart = 1;
    faStart = 1;

    % fill the sorted map with the hit trials, give the hit trials its own
    % array
    for i = 1:hitSize
        if (b{1}.hitTrialNums(i) > delTrial)
            hmhit = zeros(length(b{1}.hitTrialNums(i:hitSize)),100); % General hit array
            hmhitLick = NaN(length(b{1}.hitTrialNums(i:hitSize)),size(hm,2));
            hitNum = zeros(length(b{1}.hitTrialNums(i:hitSize)));
            hmhitpOnset = hmhit; % hit array normalized by pole onset
            hmhitpDown = hmhit; % hit array normalized by pole down
            hmhitFL = hmhit; % hit array normalized by first lick
            hmhitFAL = hmhit; % hit array normalized by first answer lick
            hmhitLL = hmhitLick; % hit array normalized by last lick
            hitDel = i;
            hitNum(1:(hitSize-hitDel+1)) = b{1}.hitTrialNums(i:hitSize)-delTrial;
            break;
        end
    end

    %Miss Trials
    missSize = length(b{1}.missTrialNums);
    missCount = 0;
    if ~(isempty(b{1}.missTrialNums))
        while (b{1}.missTrialNums(end-missCount)-delTrial >= length(start_idx)) % remove last trial if needed
            missSize = missSize-1;
            missCount = missCount+1;
        end
    end
    if missSize < 10
        missSize = 0;
    end
    % Fill sorted array with miss trials, generate miss trials array
    for i = 1:missSize
        if (b{1}.missTrialNums(i) > delTrial)
            hmmiss = zeros(length(b{1}.missTrialNums(i:missSize)),100);
            missNum = zeros(length(b{1}.missTrialNums(i:missSize)));
            hmmisspOnset = hmmiss;
            hmmisspDown = hmmiss;
            hmmissFL = hmmiss;
            missDel = i;
            missNum(1:(missSize-missDel+1)) = b{1}.missTrialNums(i:missSize)-delTrial;
            break;
        end
    end

    %Correct Rejection Trials
    CRSize = length(b{1}.correctRejectionTrialNums);
    CRCount = 0;
    if ~(isempty(b{1}.correctRejectionTrialNums))
        while (b{1}.correctRejectionTrialNums(end-CRCount)-delTrial >= length(start_idx))
            CRSize = CRSize-1;
            CRCount = CRCount+1;
        end
    end
    if CRSize < 10
        CRSize = 0;
    end
    % Fill sorted array with CR trials, generate CR trials array
    for i = 1:CRSize
        if (b{1}.correctRejectionTrialNums(i) > delTrial)
            hmCR = zeros(length(b{1}.correctRejectionTrialNums(i:CRSize)),100);
            crNum = zeros(length(b{1}.correctRejectionTrialNums(i:CRSize)));
            hmCRpOnset = hmCR;
            hmCRpDown = hmCR;
            hmCRFL =hmCR;
            crDel = i;
            crNum(1:(CRSize-crDel+1)) = b{1}.correctRejectionTrialNums(i:CRSize)-delTrial;
            break;
        end
    end

    %False Alarm Trials
    FASize = length(b{1}.falseAlarmTrialNums);
    FACount = 0;
    if ~(isempty(b{1}.falseAlarmTrialNums))
        while (b{1}.falseAlarmTrialNums(end-FACount)-delTrial >= length(start_idx))
            FASize = FASize-1;
            FACount = FACount+1;
        end
    end
    if FASize < 10
        FASize = 0;
    end
    % Fill sorted array with FA trials, generate FA trials array
    for i = 1:FASize
        if (b{1}.falseAlarmTrialNums(i) > delTrial)
            hmFA = NaN(length(b{1}.falseAlarmTrialNums(i:FASize)),size(hm,2));
            faNum = zeros(length(b{1}.falseAlarmTrialNums(i:FASize)));
            hmFApOnset = hmFA(:,1:100);
            hmFApDown = hmFA;
            hmFAFL = hmFA(:,1:100); 
            hmFAFAL = hmFA(:,1:100);
            hmFALL = hmFA;
            faDel = i;
            faNum(1:(FASize-faDel+1)) = b{1}.falseAlarmTrialNums(i:FASize)-delTrial;
            break;
        end
    end

    %% find poleonset time
    %15.44 is the frame rate
    %pull out the poleonset time from behavior file

    for i=1:length(b{1}.pinDescentOnsetTimes)-1
        poleonsetog(i)=b{1}.pinDescentOnsetTimes(i+1);
    end
    poleonset=poleonsetog(delTrial:end-1);
    poleonset=poleonset*15.44;
    for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
        polehit(i)=poleonset(b{1}.hitTrialNums(i+hitDel-1)-delTrial);
    end
    for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize))
        poleFA(i)=poleonset(b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial);
    end
    for i=1:length(b{1}.missTrialNums(missDel:missSize))
        polemiss(i)=poleonset(b{1}.missTrialNums(i+missDel-1)-delTrial);
    end
    for i=1:length(b{1}.correctRejectionTrialNums(crDel:CRSize))
        poleCR(i)=poleonset(b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial);
    end

    onset = mean(poleonset);
    frameRate = 15.44;
    secResolution = 1;

    xs = onset-frameRate:(frameRate/secResolution):size(hmhit,2); % set x axis to match
    %% find the pole down time
    for i=1:length(b{1}.pinAscentOnsetTimes)-1
        poledownog(i)=b{1}.pinAscentOnsetTimes(i+1);
    end
    poledown=poledownog(delTrial:end-1);
    poledown=poledown*15.44;
    for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
        polehitd(i)=poledown(b{1}.hitTrialNums(i+hitDel-1)-delTrial);
    end
    for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize))
        poleFAd(i)=poledown(b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial);
    end
    for i=1:length(b{1}.missTrialNums(missDel:missSize))
        polemissd(i)=poledown(b{1}.missTrialNums(i+missDel-1)-delTrial);
    end
    for i=1:length(b{1}.correctRejectionTrialNums(crDel:CRSize))
        poleCRd(i)=poledown(b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial);
    end

    %% find and rasterplot first lick on the heatmap
    for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
        firstlickhit(i)=b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes(1);
    end
    firstlickhit=firstlickhit*15.44;

    for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize))
        firstlickFA(i)=b{1}.trials{b{1}.falseAlarmTrialNums(i+faDel-1)}.beamBreakTimes(1);
    end
    
    if FASize ~= 0
        firstlickFA=firstlickFA*15.44;
    end
    
    for i=1:length(b{1}.missTrialNums(missDel:missSize))
        if isempty(b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes)==1
            firstlickMiss(i)=0;
        else
            if b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes(1)>length(tmp2{b{1}.missTrialNums(i+missDel-1)-delTrial})/15.44
              firstlickMiss(i)=0;
            else
            firstlickMiss(i)=b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes(1);
            end
        end
    end
    if missSize ~= 0
        firstlickMiss=firstlickMiss*15.44;
    end
    
    for i=1:length(b{1}.correctRejectionTrialNums(crDel:CRSize))
        if isempty(b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes)==1
            firstlickCR(i)=0;
        else
            if b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes(1)>length(tmp2{b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial})/15.44
            firstlickCR(i)=0;
            else
            firstlickCR(i)=b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes(1);
            end
        end
    end
    
    if CRSize ~= 0
        firstlickCR=firstlickCR*15.44;
    end




    %% find and rasterplot the first answer lick time
    for i=1:length(b{1}.hitTrialNums(hitDel:hitSize)) %Loop through trials, find the first time the beam is broken(lick) that occurs after the answer period begins
        FALTHitVal = min(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes > b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.answerPeriodTime(1)));
        if isempty(FALTHitVal)
            FALTHit(i) = 0;
        elseif FALTHitVal > b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.answerPeriodTime(2)
            FALTHitVal = min(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes > b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.answerPeriodTime(1)-0.1)); % Correction by one millisecond if needed
            FALTHit(i) = FALTHitVal;
        else
            FALTHit(i) = FALTHitVal;
        end
    end
    FALTHit= FALTHit*15.44; % multiply by frame rate
    for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize)) %Loop, find first break during false alarm trials
        FALTFAVal=min(b{1}.trials{b{1}.falseAlarmTrialNums(i+faDel-1)}.beamBreakTimes(b{1}.trials{b{1}.falseAlarmTrialNums(i+faDel-1)}.beamBreakTimes >b{1}.trials{b{1}.falseAlarmTrialNums(i+faDel-1)}.answerPeriodTime(1)));
        if isempty(FALTFAVal)
            FALTFA(i) = 0;
        else
            FALTFA(i) = FALTFAVal;
        end
    end
    
    if FASize ~= 0
        FALTFA=FALTFA*15.44;
    end
    
    %% find last lick
    %     for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
    %         lastlickhit(i)=b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes(end);
    %     end

        for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
            lastlickhit(i)=max(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes < (length(tmp2{b{1}.hitTrialNums(i+hitDel-1)-delTrial})/15.44)));
        end   
        lastlickhit=lastlickhit*15.44;

        for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize))
            lastlickFA(i)=max(b{1}.trials{b{1}.falseAlarmTrialNums(i+faDel-1)}.beamBreakTimes(b{1}.trials{b{1}.falseAlarmTrialNums(i+faDel-1)}.beamBreakTimes < (length(tmp2{b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial})/15.44)));
        end   
        
        if FASize ~= 0
            lastlickFA=lastlickFA*15.44;
        end
        
        for i=1:length(b{1}.missTrialNums(missDel:missSize))
    %         if isempty(b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes)==1
            if isempty(max(b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes(b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes < (length(tmp2{b{1}.missTrialNums(i+missDel-1)-delTrial})/15.44))))==1 
                lastlickMiss(i)=0;
            else
                lastlickMiss(i)=max(b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes(b{1}.trials{b{1}.missTrialNums(i+missDel-1)}.beamBreakTimes < (length(tmp2{b{1}.missTrialNums(i+missDel-1)-delTrial})/15.44)));
            end
        end
        if missSize ~= 0
            lastlickMiss=lastlickMiss*15.44;
        end
        
        for i=1:length(b{1}.correctRejectionTrialNums(crDel:CRSize))
    %         if isempty(b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes)==1
            if isempty(max(b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes(b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes < (length(tmp2{b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial})/15.44))))==1
                lastlickCR(i)=0;
            else
                lastlickCR(i)=max(b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes(b{1}.trials{b{1}.correctRejectionTrialNums(i+crDel-1)}.beamBreakTimes < (length(tmp2{b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial})/15.44)));
            end
        end
        if CRSize ~= 0
            lastlickCR=lastlickCR*15.44;
        end

    %% Add hm based on 200 ms. before 
    timeFrame = 0.5 * frameRate;

    %% loop EVERYTHING to be with normalization
    for loopcount = 1:3
        if loopcount == 1
            hmuse = hm2;
            hmuseFA = hm;
            hmuseHit = hm;
            meanFolder = timeMeanFolder;
            moveFolder = moveFolder3;
        elseif loopcount == 2
            hmuse = hm4;
            hmuseFA = hm3;
            hmuseHit = hm3;
            meanFolder = normMeanFolder;
            moveFolder = moveFolder2;
        else
            hmuse = hm2;
            hmuseFA = hm;
            hmuseHit = hm;
            meanFolder = baseMeanFolder;
            moveFolder = moveFolder1;
        end
        %% plot the sorted heatmap by trial type
        %Automated Sorting by trial type with average responses ready to be plotted

        % initialize counters
        counter1 = 1;
        counter2 = 1;
        % fill the sorted map with the hit trials, give the hit trials its own
        % array



         for i = hitDel:hitSize
            if (b{1}.hitTrialNums(i) > delTrial)
                for (k = i:hitSize)
                    hmbsort(counter1,:) = hmuse((b{1}.hitTrialNums(k)-delTrial),:);
                    hmhit(counter2,:) = hmuse((b{1}.hitTrialNums(k)-delTrial),:);
                    hmhitLick(counter2,:) = hmuseHit((b{1}.hitTrialNums(k)-delTrial),:);

                    % If normalization by time, normalize by all prior
                    % measures, else, all will be original values
                    if loopcount == 1
                        meanuse = mean(tmp2{hitNum(k-hitDel+1)}((floor(polehit(k-hitDel+1)-timeFrame)):floor(polehit(k-hitDel+1))));
                        if (isnan(meanuse))
                            hmhitpOnset(counter2,:) = hmhit(counter2,:);
                        else
                            hmhitpOnset(counter2,:) = (hmhit(counter2,:)+1)/meanuse*(tmp3(hitNum(k-hitDel+1)))-1;
                        end

                        meanuse = mean(tmp2{hitNum(k-hitDel+1)}((floor(polehitd(k-hitDel+1)-timeFrame)):floor(polehitd(k-hitDel+1)))); 
                        if (isnan(meanuse))
                            hmhitpDown(counter2,:) = hmhit(counter2,:);
                        else
                            hmhitpDown(counter2,:) = (hmhit(counter2,:)+1)/meanuse*(tmp3(hitNum(k-hitDel+1)))-1;
                        end

                        if (floor(firstlickhit(k-hitDel+1)-timeFrame)) <= 0
                            timeStart = 1;
                        else
                            timeStart = floor(firstlickhit(k-hitDel+1)-timeFrame);
                        end

                        meanuse = mean(tmp2{hitNum(k-hitDel+1)}(timeStart:floor(firstlickhit(k-hitDel+1))));
                        if (isnan(meanuse))
                            hmhitFL(counter2,:) = hmhit(counter2,:);
                        else
                            hmhitFL(counter2,:) = (hmhit(counter2,:)+1)/meanuse*(tmp3(hitNum(k-hitDel+1)))-1;
                        end

                        if (FALTHit(k-hitDel+1)-timeFrame) <= 0
                            timeStart = 1;
                        else
                            timeStart = floor(FALTHit(k-hitDel+1)-timeFrame);
                        end

                        meanuse = mean(tmp2{hitNum(k-hitDel+1)}(timeStart:floor(FALTHit(k-hitDel+1))));
                        if (isnan(meanuse))
                            hmhitFAL(counter2,:) = hmhit(counter2,:);
                        else
                            hmhitFAL(counter2,:) = (hmhit(counter2,:)+1)/meanuse*(tmp3(hitNum(k-hitDel+1)))-1;
                        end

                        meanuse = mean(tmp2{hitNum(k-hitDel+1)}((floor(lastlickhit(k-hitDel+1)-timeFrame)):floor(lastlickhit(k-hitDel+1))));
                        if (isnan(meanuse))
                            hmhitLL(counter2,:) = hmhit(counter2,:);
                        else
                            hmhitLL(counter2,:) = (hmhitLick(counter2,:)+1)/meanuse*(tmp3(hitNum(k-hitDel+1)))-1;
                        end

                    else
                        hmhitpOnset = hmhit;
                        hmhitpDown = hmhit;
                        hmhitFL = hmhit;
                        hmhitFAL = hmhit;
                        hmhitLL = hmhitLick;
                    end

                    counter1 = counter1+1;
                    counter2 = counter2+1;
                end
                break;
            end
        end
        counter1 = counter1 + 50; % space between the hit and miss trials

        %Miss Trials
        missStart = counter1;
        counter2 = 1;

        % Fill sorted array with miss trials, generate miss trials array

    %     for i = 1:missSize
        for i = missDel:missSize
            if (b{1}.missTrialNums(i) > delTrial)
                for (k = i:missSize)
                    hmbsort(counter1,:) = hmuse((b{1}.missTrialNums(k)-delTrial),:);
                    hmmiss(counter2,:) = hmuse((b{1}.missTrialNums(k)-delTrial),:);
                    % If normalization by time, normalize by all prior
                    % measures, else, all will be original values
                    if loopcount == 1
                        meanuse = mean(tmp2{missNum(k-missDel+1)}((floor(polemiss(k-missDel+1)-timeFrame)):floor(polemiss(k-missDel+1))));
                        if (isnan(meanuse))
                            hmmisspOnset(counter2,:) = hmmiss(counter2,:);
                        else
                            hmmisspOnset(counter2,:) = (hmmiss(counter2,:)+1)/meanuse*(tmp3(missNum(k-missDel+1)))-1;
                        end                  

                        meanuse = mean(tmp2{missNum(k-missDel+1)}((floor(polemissd(k-missDel+1)-timeFrame)):floor(polemissd(k-missDel+1))));
                        if (isnan(meanuse))
                            hmmisspDown(counter2,:) = hmmiss(counter2,:);
                        else
                            hmmisspDown(counter2,:) = (hmmiss(counter2,:)+1)/meanuse*(tmp3(missNum(k-missDel+1)))-1;
                        end                                    

                        if (floor(firstlickMiss(k-missDel+1)-timeFrame)) <= 0
                            timeStart = 1;
                        else
                            timeStart = floor(firstlickMiss(k-missDel+1)-timeFrame);
                        end    

                        if (firstlickMiss(k-missDel+1)) ~= 0
                            meanuse = mean(tmp2{missNum(k-missDel+1)}(timeStart:floor(firstlickMiss(k-missDel+1))));
                        else
                            meanuse = nan;
                        end

                        if (isnan(meanuse))
                            hmmissFL(counter2,:) = hmmiss(counter2,:);
                        else
                            hmmissFL(counter2,:) = (hmmiss(counter2,:)+1)/meanuse*(tmp3(missNum(k-missDel+1)))-1;
                        end                        

                    else
                        hmmisspOnset = hmmiss;
                        hmmisspDown = hmmiss;
                        hmmissFL = hmmiss;
                    end
                    counter1 = counter1+1;
                    counter2 = counter2+1;
                end
                break;
            end
        end
        counter1 = counter1 + 50;

        %Correct Rejection Trials
        crStart = counter1;
        counter2 = 1;

        % Fill sorted array with CR trials, generate CR trials array
    %     for i = 1:CRSize
        for i = crDel:CRSize
            if (b{1}.correctRejectionTrialNums(i) > delTrial)
                for (k = i:CRSize)
                    hmbsort(counter1,:) = hmuse((b{1}.correctRejectionTrialNums(k)-delTrial),:);
                    hmCR(counter2,:) = hmuse((b{1}.correctRejectionTrialNums(k)-delTrial),:);
                    if loopcount == 1
                        meanuse = mean(tmp2{crNum(k-crDel+1)}((floor(poleCR(k-crDel+1)-timeFrame)):floor(poleCR(k-crDel+1))));
                        if (isnan(meanuse))
                            hmCRpOnset(counter2,:) = hmCR(counter2,:);
                        else
                            hmCRpOnset(counter2,:) = (hmCR(counter2,:)+1)/meanuse*(tmp3(crNum(k-crDel+1)))-1;
                        end                         

                        meanuse = mean(tmp2{crNum(k-crDel+1)}((floor(poleCRd(k-crDel+1)-timeFrame)):floor(poleCRd(k-crDel+1))));
                        if (isnan(meanuse))
                            hmCRpDown(counter2,:) = hmCR(counter2,:);
                        else
                            hmCRpDown(counter2,:) = (hmCR(counter2,:)+1)/meanuse*(tmp3(crNum(k-crDel+1)))-1;
                        end                                     

                        if (floor(firstlickCR(k-crDel+1)-timeFrame)) <= 0
                            timeStart = 1;
                        else
                            timeStart = floor(firstlickCR(k-crDel+1)-timeFrame);
                        end        

                        if (firstlickCR(k-crDel+1)) ~= 0
                            meanuse = mean(tmp2{crNum(k-crDel+1)}(timeStart:floor(firstlickCR(k-crDel+1))));
                        else
                            meanuse = nan;
                        end

                        if (isnan(meanuse))
                            hmCRFL(counter2,:) = hmCR(counter2,:);
                        else
                            hmCRFL(counter2,:) = (hmCR(counter2,:)+1)/meanuse*(tmp3(crNum(k-crDel+1)))-1;
                        end                   

                    else
                        hmCRpOnset = hmCR;
                        hmCRpDown = hmCR;
                        hmCRFL= hmCR;
                    end
                    counter1 = counter1+1;
                    counter2 = counter2+1;
                end
                break;
            end
        end
            counter1 = counter1 + 50;

        %False Alarm Trials
        faStart = counter1;
        counter2 = 1;

        % Fill sorted array with FA trials, generate FA trials array
    %     for i = 1:FASize
        for i = faDel:FASize
            if (b{1}.falseAlarmTrialNums(i) > delTrial)
                for (k = i:FASize)
                    hmbsort(counter1,:) = hmuse((b{1}.falseAlarmTrialNums(k)-delTrial),:);
                    hmFA(counter2,:) = hmuseFA((b{1}.falseAlarmTrialNums(k)-delTrial),:);
                    if loopcount == 1
                        meanuse = mean(tmp2{faNum(k-faDel+1)}((floor(poleFA(k-faDel+1)-timeFrame)):floor(poleFA(k-faDel+1))));
                        if (isnan(meanuse))
                            hmFApOnset(counter2,1:100) = hmFA(counter2,1:100);
                        else
                            hmFApOnset(counter2,1:100) = (hmFA(counter2,1:100)+1)/meanuse*(tmp3(faNum(k-faDel+1)))-1;
                        end 

                        meanuse = mean(tmp2{faNum(k-faDel+1)}((floor(poleFAd(k-faDel+1)-timeFrame)):floor(poleFAd(k-faDel+1))));
                        if (isnan(meanuse))
                            hmFApDown(counter2,:) = hmFA(counter2,:);
                        else
                            hmFApDown(counter2,:) = (hmFA(counter2,:)+1)/meanuse*(tmp3(faNum(k-faDel+1)))-1;
                        end 

                        if (floor(firstlickFA(k-faDel+1)-timeFrame)) <= 0
                            timeStart = 1;
                        else
                            timeStart = floor(firstlickFA(k-faDel+1)-timeFrame);            
                        end

                        meanuse = mean(tmp2{faNum(k-faDel+1)}(timeStart:floor(firstlickFA(k-faDel+1))));
                        if (isnan(meanuse))
                            hmFAFL(counter2,1:100) = hmFA(counter2,1:100);
                        else
                            hmFAFL(counter2,1:100) = (hmFA(counter2,1:100)+1)/meanuse*(tmp3(faNum(k-faDel+1)))-1;
                        end

                        if (floor(FALTFA(k-faDel+1)-timeFrame)) <= 0
                            timeStart = 1;
                        else
                            timeStart = floor(FALTFA(k-faDel+1)-timeFrame);
                        end

                        if FALTFA(k-faDel+1) > length(tmp2{faNum(k-faDel+1)})
                            meanuse = mean(tmp2{faNum(k-faDel+1)}(timeStart:length(tmp2{faNum(k-faDel+1)})));
                        else
                            meanuse = mean(tmp2{faNum(k-faDel+1)}(timeStart:floor(FALTFA(k-faDel+1))));
                        end

                        if (isnan(meanuse))
                            hmFAFAL(counter2,1:100) = hmFA(counter2,1:100);
                        else
                            hmFAFAL(counter2,1:100) = (hmFA(counter2,1:100)+1)/meanuse*(tmp3(faNum(k-faDel+1)))-1;
                        end

                        meanuse = mean(tmp2{faNum(k-faDel+1)}((floor(lastlickFA(k-faDel+1)-timeFrame)):floor(lastlickFA(k-faDel+1))));
                        if (isnan(meanuse))
                            hmFALL(counter2,:) = hmFA(counter2,:);
                        else
                            hmFALL(counter2,:) = (hmFA(counter2,:)+1)/meanuse*(tmp3(faNum(k-faDel+1)))-1;
                        end 

                    else
                        hmFApOnset = hmFA(:,1:100);
                        hmFApDown = hmFA;
                        hmFAFL = hmFA(:,1:100);
                        hmFAFAL = hmFA(:,1:100);
                        hmFALL = hmFA;
                    end
                    counter1 = counter1+1;
                    counter2 = counter2+1;
                end
                break;
            end
        end

        counter1 = counter1 + 50;


        %% Visualize the sorted heat map

        %normalization by choosing the first 16 frames as the baseline,
        %subtracting the baseline from the original
        hitbase=mean(nanmean(hmhitpOnset(:,1:16)));
        if missSize ~= 0
            missbase=mean(nanmean(hmmisspOnset(:,1:16)));
        end
        if CRSize ~= 0
            CRbase=mean(nanmean(hmCRpOnset(:,1:16)));
        end
        if FASize ~= 0
            FAbase=mean(nanmean(hmFApOnset(:,1:16)));
        end
        normhit=(nanmean(hmhitpOnset)-hitbase);
        if missSize ~= 0
            normmiss=(nanmean(hmmisspOnset)-missbase);
        end
        if CRSize ~= 0
            normCR=(nanmean(hmCRpOnset)-CRbase);
        end
        if FASize ~= 0
            normFA=(nanmean(hmFApOnset)-FAbase);
        end
 
        %% rasterplot all the lick on the heatmap for different trial type
        %plot all the lick for hit trials
        for i=1:length(b{1}.hitTrialNums(hitDel:hitSize)) % get all lick indices
            hitlicknum(i)=numel(b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes);
        end

        hitlick=zeros(length(b{1}.hitTrialNums(hitDel:hitSize)),max(hitlicknum));
        for i=1:length(b{1}.hitTrialNums(hitDel:end)) % get all licks
            lickT = b{1}.trials{b{1}.hitTrialNums(i+hitDel-1)}.beamBreakTimes;
            hitlick(i,1:length(lickT))=(lickT)*15.44;
        end


        %% hit trials aligned to first answer lick
        rndFALTHit=ceil(max(FALTHit));
        lickalign=nan(length(b{1}.hitTrialNums(hitDel:hitSize)),100+rndFALTHit); % Changed here end-1
        for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
            tSize = length(hmuse((b{1}.hitTrialNums(i+hitDel-1)-delTrial),:))-1;
            shift = ceil(max(FALTHit) - FALTHit(i)) + 1;
            lickalign(i,shift:(shift+tSize))=hmhitFAL(i,:); % align trials to the furthest up answer lick across trials, align all answer licks
        end
        for i=1:length(FALTHit)
            maxFALTHit(i)=max(FALTHit);
        end

        % Plot heat map with first answer lick
        
        FALTHitAlign = max(maxFALTHit) * ones(1,2);
        if (size(lickalign,1) == 1)
            responsetrace = lickalign;
        else
            responsetrace = nanmean(lickalign);
        end

        % plot average intensity trace
        xs3 = max(maxFALTHit)-frameRate:(frameRate/secResolution):size(lickalign,2);
        currxlim = [max(maxFALTHit)-1.5*frameRate, max(maxFALTHit)+4.5*frameRate];
        if loopcount == 1
            currxlim = [max(maxFALTHit)-timeFrame, max(maxFALTHit)+2.05*frameRate];
        end
        responsetraceview = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace(abs(ceil(currxlim(1))):ceil(currxlim(2))));
        responsetraceview(1:abs(floor(currxlim(1)))) = nan;
        
        hitAlignFAL = responsetraceview(abs(ceil(currxlim(1))):ceil(currxlim(2)));
        %aligned heatmap with average response Lily0906
 
 
        % use 200ms before first reward lick as baseline


        %% FA aligned to first answer lick
        if FASize ~= 0
            rndFALTFA=ceil(max(FALTFA));
            lickalign3=nan(length(b{1}.falseAlarmTrialNums(faDel:FASize)),100+rndFALTFA);
            for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize)) % align answer licks
                tSize3 = length(hmuse((b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial),:))-1;
                shift3 = ceil(max(FALTFA) - FALTFA(i)) + 1;
                lickalign3(i,shift3:(shift3+tSize3))=hmFAFAL(i,1:100);
            end
            for i=1:length(FALTFA)
                maxFALTFA(i)=max(FALTFA);
            end

            % plot heat map

            FALTFAAlign = max(maxFALTFA) * ones(1,2);
            if (size(lickalign3,1) == 1)
                responsetrace3 = lickalign3;
            else
                responsetrace3= nanmean(lickalign3);
            end
            xs = max(maxFALTFA)-frameRate:(frameRate/secResolution):size(lickalign3,2);
            currxlim = [max(maxFALTFA)-1.5*frameRate, max(maxFALTFA)+4.5*frameRate];
            if loopcount == 1
                currxlim = [max(maxFALTFA)-timeFrame,max(maxFALTFA)+2.05*frameRate];
            end
            % plot average intensity trace
            responsetrace3view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace3(abs(ceil(currxlim(1))):ceil(currxlim(2))));
            responsetrace3view(1:abs(floor(currxlim(1)))) = nan;
            FAAlignFAL = responsetrace3view(abs(ceil(currxlim(1))):ceil(currxlim(2)));
        end
        %% hit trials aligned to first lick
        rndfirstlickhit=ceil(max(firstlickhit));
        lickalign2=nan(length(b{1}.hitTrialNums(hitDel:hitSize)),100+rndfirstlickhit);
        for i=1:length(b{1}.hitTrialNums(hitDel:hitSize)) % align licks
            tSize2 = length(hmuse((b{1}.hitTrialNums(i+hitDel-1)-delTrial),:))-1;
            shift2 = ceil(max(firstlickhit) - firstlickhit(i)) + 1;
            lickalign2(i,shift2:(shift2+tSize2))=hmhitFL(i,:);
        end
        for i=1:length(firstlickhit)
            maxfirstlickhit(i)=max(firstlickhit);
        end

        % plot heat map
  
        if (size(lickalign2,1) == 1)
            responsetrace2 = lickalign2;
        else
            responsetrace2= nanmean(lickalign2);
        end

        % plot average intensity trace
        xs = max(maxfirstlickhit)-frameRate:(frameRate/secResolution):size(lickalign2,2);
        currxlim = [max(maxfirstlickhit)-1.5*frameRate, max(maxfirstlickhit)+4.5*frameRate];
        if loopcount == 1
            currxlim = [max(maxfirstlickhit)-timeFrame, max(maxfirstlickhit)+2.05*frameRate];
        end
        responsetrace2view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace2(abs(ceil(currxlim(1))):ceil(currxlim(2))));
        responsetrace2view(1:abs(floor(currxlim(1)))) = nan;
        hitAlignFL = responsetrace2view(abs(ceil(currxlim(1))):ceil(currxlim(2)));

        FLHitAlign = max(maxfirstlickhit) * ones(1,2);
        
        %% FA trials aligned to first lick
        if FASize ~= 0
            rndfirstlickFA=ceil(max(firstlickFA));
            lickalign4=nan(length(b{1}.falseAlarmTrialNums(faDel:FASize)),100+rndfirstlickFA);
            for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize)) % align licks
                tSize4 = length(hmuse((b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial),:))-1;
                shift4 = ceil(max(firstlickFA) - firstlickFA(i)) + 1;
                lickalign4(i,shift4:(shift4+tSize4))=hmFAFL(i,1:100);
            end
            for i=1:length(firstlickFA)
                maxfirstlickFA(i)=max(firstlickFA);
            end

            % plot heat map

            if (size(lickalign4,1) == 1)
                responsetrace4 = lickalign4;
            else
                responsetrace4= nanmean(lickalign4);
            end

            % plot average intensity trace

            xs = max(maxfirstlickFA)-frameRate:(frameRate/secResolution):size(lickalign4,2);
            currxlim = [max(maxfirstlickFA)-1.5*frameRate, max(maxfirstlickFA)+4.5*frameRate];
            if loopcount == 1
                currxlim = [max(maxfirstlickFA)-timeFrame, max(maxfirstlickFA)+2.05*frameRate];
            end
            responsetrace4view = horzcat(zeros(1,abs(floor(xs(1)))),responsetrace4(abs(ceil(xs(1))):ceil(currxlim(2))));
            responsetrace4view(1:abs(floor(currxlim(1)))) = nan;
            FAAlignFL = responsetrace4view(abs(ceil(currxlim(1))):ceil(currxlim(2)));
            FLFAAlign = max(maxfirstlickFA) * ones(1,2);
        end

        

         %% Miss trial aligned to first lick
        if missSize ~= 0
            rndfirstlickMiss=ceil(max(firstlickMiss));
            lickalign12=nan(length(b{1}.missTrialNums(missDel:missSize)),100+rndfirstlickMiss);
                for i=1:length(b{1}.missTrialNums(missDel:missSize)) % align licks
                    tSize12 = length(hmuse((b{1}.missTrialNums(i+missDel-1)-delTrial),:))-1;
                    shift12 = ceil(max(firstlickMiss) - firstlickMiss(i)) + 1;
                    lickalign12(i,shift12:(shift12+tSize12))=hmmissFL(i,:);
                end
                for i=1:length(firstlickMiss)
                    maxfirstlickMiss(i)=max(firstlickMiss);
                end
            firstlickMiss2 = firstlickMiss(~(firstlickMiss == 0));
            lickalign12true=lickalign12(~(firstlickMiss == 0), :);  
            maxfirstlickMiss2=maxfirstlickMiss(1:length(firstlickMiss2));     
            if ~(isempty(maxfirstlickMiss2))        
                % plot heat map

                if size(lickalign12true,1) == 1
                    responsetrace12 = lickalign12true;
                else
                    responsetrace12 = nanmean(lickalign12true);
                end
              % plot average intensity trace

                xs12 = max(maxfirstlickMiss)-frameRate:(frameRate/secResolution):size(lickalign12,2);
                currxlim = [max(maxfirstlickMiss)-1.5*frameRate, max(maxfirstlickMiss)+4.5*frameRate];
                if loopcount == 1
                    currxlim = [max(maxfirstlickMiss)-timeFrame, max(maxfirstlickMiss)+2.05*frameRate];
                end
                if (currxlim(1) > 0)
                    responsetrace12view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace12(abs(ceil(currxlim(1))):ceil(currxlim(2))));
                    responsetrace12view(1:abs(floor(currxlim(1)))) = nan;
                    missAlignFL = responsetrace12view(abs(ceil(currxlim(1))):ceil(currxlim(2)));
                else
                    responsetrace12view = horzcat(zeros(1,abs(floor(xs12(1)))),responsetrace12(abs(ceil(xs12(1))):ceil(currxlim(2))));
                    responsetrace12view(1:abs(floor(xs12(1)))) = nan;                    
                    missAlignFL = responsetrace12view(abs(ceil(xs12(1))):ceil(currxlim(2)));
                end


                FLMissAlign = max(maxfirstlickMiss2) * ones(1,2);
            end
        end
         %% CR trials aligned to first lick
         if CRSize ~= 0
            rndfirstlickCR=ceil(max(firstlickCR));
            lickalign13=nan(length(b{1}.correctRejectionTrialNums(crDel:CRSize)),100+rndfirstlickCR);
                for i=1:length(b{1}.correctRejectionTrialNums(crDel:CRSize)) % align licks
                    tSize13 = length(hmuse((b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial),:))-1;
                    shift13 = ceil(max(firstlickCR) - firstlickCR(i)) + 1;
                    lickalign13(i,shift13:(shift13+tSize13))=hmCRFL(i,:);
                end
                for i=1:length(firstlickCR)
                    maxfirstlickCR(i)=max(firstlickCR);
                end
            firstlickCR2 = firstlickCR(~(firstlickCR == 0));
            lickalign13true=lickalign13(~(firstlickCR == 0), :);  
            maxfirstlickCR2=maxfirstlickCR(1:length(firstlickCR2));     
            if ~(isempty(maxfirstlickCR2))        
                % plot heat map

                if (size(lickalign13true,1) == 1)
                    responsetrace13 = lickalign13true;
                else
                    responsetrace13= nanmean(lickalign13true);
                end
                % plot average intensity trace
                xs13 = max(maxfirstlickCR)-frameRate:(frameRate/secResolution):size(lickalign13,2);
                currxlim = [max(maxfirstlickCR)-1.5*frameRate, max(maxfirstlickCR)+4.5*frameRate];
                if loopcount == 1
                    currxlim = [max(maxfirstlickCR)-timeFrame, max(maxfirstlickCR)+2.05*frameRate];
                end
                if (currxlim(1) > 0)
                    responsetrace13view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace13(abs(ceil(currxlim(1))):ceil(currxlim(2))));
                    responsetrace13view(1:abs(floor(currxlim(1)))) = nan;
                    CRAlignFL = responsetrace13view(abs(ceil(currxlim(1))):ceil(currxlim(2)));
                else
                    responsetrace13view = horzcat(zeros(1,abs(floor(xs13(1)))),responsetrace13(abs(ceil(xs13(1))):ceil(currxlim(2))));
                    responsetrace13view(1:abs(floor(xs13(1)))) = nan;
                    CRAlignFL = responsetrace13view(abs(ceil(xs13(1))):ceil(currxlim(2)));
                end
                FLCRAlign = max(maxfirstlickCR2) * ones(1,2);
            end
         end


        %% aligned to the first lick using the first lick is not first reward lick

        firstlickhit2 = firstlickhit(~(firstlickhit == FALTHit));

        lickalign2true=lickalign2(~(firstlickhit == FALTHit),:);

        maxfirstlickhit2=maxfirstlickhit(1:length(firstlickhit2));

        % plot heat map
        
        responsetrace2true= nanmean(lickalign2true);

        % plot average intensity trace
        xs = max(maxfirstlickhit2)-frameRate:(frameRate/secResolution):size(lickalign2true,2);
        currxlim = [max(maxfirstlickhit2)-1.5*frameRate, max(maxfirstlickhit2)+4.5*frameRate];
        if loopcount == 1
            currxlim = [max(maxfirstlickhit2)-timeFrame, max(maxfirstlickhit2)+2.05*frameRate];
        end
        if (currxlim(1) > 0)
            responsetrace2viewtrue = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace2true(abs(ceil(currxlim(1))):ceil(currxlim(2))));
            responsetrace2viewtrue(1:abs(floor(currxlim(1)))) = nan;
            hitAlignFLtrue = responsetrace2viewtrue(abs(ceil(currxlim(1))):ceil(currxlim(2)));
        else 
            responsetrace2viewtrue = horzcat(zeros(1,abs(floor(xs(1)))),responsetrace2true(abs(ceil(xs(1))):ceil(currxlim(2))));
            responsetrace2viewtrue(1:abs(floor(xs(1)))) = nan;
            hitAlignFLtrue = responsetrace2viewtrue(abs(ceil(xs(1))):ceil(currxlim(2)));
        end       
      
        FLHitAlign2 = max(maxfirstlickhit2) * ones(1,2);
        

        %% hit trials align to pole down
        rndpoledownhit=ceil(max(polehitd));
        pdalign=nan(length(b{1}.hitTrialNums(hitDel:hitSize)),100+rndpoledownhit);
        for i=1:length(b{1}.hitTrialNums(hitDel:hitSize))
            tSize1 = length(hmuse((b{1}.hitTrialNums(i+hitDel-1)-delTrial),:))-1;
            shift5 = ceil(max(polehitd) - polehitd(i)) + 1;
            pdalign(i,shift5:(shift5+tSize1))=hmhitpDown(i,:); % align to pole down
        end
        for i=1:length(polehitd)
            maxpolehitd(i)=max(polehitd);
        end

       
        responsetrace5 = nanmean(pdalign);

        % plot average response trace
        xs = max(maxpolehitd)-2*frameRate:(frameRate/secResolution):size(pdalign,2);
        currxlim = [max(maxpolehitd)-2.5*frameRate, max(maxpolehitd)+2.5*frameRate];
        if loopcount == 1
            currxlim = [max(maxpolehitd)-timeFrame, max(maxpolehitd)+2.05*frameRate];
        end
        responsetrace5view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace5(abs(ceil(currxlim(1))):ceil(currxlim(2))));

        responsetrace5view(1:abs(floor(currxlim(1)))) = nan;
        
        hitAlignPD = responsetrace5view(abs(ceil(currxlim(1))):ceil(currxlim(2)));

        PDHitAlign = max(maxpolehitd) * ones(1,2);
        
        %% FA trials align to pole down
        if FASize ~= 0
            rndpoledownFA=ceil(max(poleFAd));
            pdalign4=nan(length(b{1}.falseAlarmTrialNums(faDel:FASize)),100+rndpoledownFA);
            for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize)) % align to pole down
                tSize4 = length(hmuseFA((b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial),:))-1;
                shift4 = ceil(max(poleFAd) - poleFAd(i)) + 1;
                pdalign4(i,shift4:(shift4+tSize4))=hmFApDown(i,:);
            end
            for i=1:length(poleFAd)
                maxpoleFAd(i)=max(poleFAd);
            end

            % plot heat map

            responsetrace8= nanmean(pdalign4);
            xs = max(maxpoleFAd)-2*frameRate:(frameRate/secResolution):size(pdalign4,2);
            currxlim = [max(maxpoleFAd)-2.5*frameRate, max(maxpoleFAd)+2.5*frameRate];
            if loopcount == 1
                currxlim = [max(maxpoleFAd)-timeFrame, max(maxpoleFAd)+2.05*frameRate];
            end
            responsetrace8view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace8(abs(ceil(currxlim(1))):ceil(currxlim(2))));

            responsetrace8view(1:abs(floor(currxlim(1)))) = nan;
            
            FAAlignPD = responsetrace8view(abs(ceil(currxlim(1))):ceil(currxlim(2)));

            PDFAAlign = max(maxpoleFAd) * ones(1,2);
        end
        
        %% miss trials align to pole down
        if missSize ~= 0
            rndpoledownmiss=ceil(max(polemissd));
            pdalign2=nan(length(b{1}.missTrialNums(missDel:missSize)),100+rndpoledownmiss);   %0918 Lily before was end
            for i=1:length(b{1}.missTrialNums(missDel:missSize))
                tSize7 = length(hmuse((b{1}.missTrialNums(i+missDel-1)-delTrial),:))-1;
                shift7 = ceil(max(polemissd) - polemissd(i)) + 1;
                pdalign2(i,shift7:(shift7+tSize7))=hmmisspDown(i,:);
            end
            for i=1:length(polemissd)
                maxpolemissd(i)=max(polemissd);
            end

            responsetrace6= nanmean(pdalign2);
            xs = max(maxpolemissd)-2*frameRate:(frameRate/secResolution):size(pdalign2,2);
            currxlim = [max(maxpolemissd)-2.5*frameRate, max(maxpolemissd)+2.5*frameRate];
            if loopcount == 1
                currxlim = [max(maxpolemissd)-timeFrame, max(maxpolemissd)+2.05*frameRate];
            end
            responsetrace6view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace6(abs(ceil(currxlim(1))):ceil(currxlim(2))));

            responsetrace6view(1:abs(floor(currxlim(1)))) = nan;
            
            missAlignPD = responsetrace6(abs(ceil(currxlim(1))):ceil(currxlim(2)));

            PDMissAlign = max(maxpolemissd) * ones(1,2);
        end
        
        %% CR trials align to pole down
        if CRSize ~= 0
            rndpoledownCR=ceil(max(poleCRd));
            pdalign3=nan(length(b{1}.correctRejectionTrialNums(crDel:CRSize)),100+rndpoledownCR);    %Lily0918 before was to the end
            for i=1:length(b{1}.correctRejectionTrialNums(crDel:CRSize))   %Lily0918 has -1 before
                tSize7 = length(hmuse((b{1}.correctRejectionTrialNums(i+crDel-1)-delTrial),:))-1;
                shift7 = ceil(max(poleCRd) - poleCRd(i)) + 1;
                pdalign3(i,shift7:(shift7+tSize7))=hmCRpDown(i,:);
            end
            for i=1:length(poleCRd)
                maxpoleCRd(i)=max(poleCRd);
            end


            responsetrace7= nanmean(pdalign3);
            xs = max(maxpoleCRd)-2*frameRate:(frameRate/secResolution):size(pdalign3,2);
            currxlim = [max(maxpoleCRd)-2.5*frameRate, max(maxpoleCRd)+2.5*frameRate];
            if loopcount == 1
                currxlim = [max(maxpoleCRd)-timeFrame, max(maxpoleCRd)+2.05*frameRate];
            end
            responsetrace7view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace7(abs(ceil(currxlim(1))):ceil(currxlim(2))));

            responsetrace7view(1:abs(floor(currxlim(1)))) = nan;
            
            CRAlignPD = responsetrace7view(abs(ceil(currxlim(1))):ceil(currxlim(2)));

            PDCRAlign = max(maxpoleCRd) * ones(1,2);
        end

        %% hit trials aligned to last lick
        rndlastlickhit=ceil(max(lastlickhit));
        lickalign10=NaN(length(b{1}.hitTrialNums(hitDel:hitSize)),100+rndlastlickhit);
        for i=1:length(b{1}.hitTrialNums(hitDel:hitSize)) % align licks
            tSize10 = length(hmuseHit((b{1}.hitTrialNums(i+hitDel-1)-delTrial),:))-1;
            shift10 = ceil(max(lastlickhit) - lastlickhit(i)) + 1;
            lickalign10(i,shift10:(shift10+tSize10))=hmhitLL(i,:);
        end
        for i=1:length(lastlickhit)
            maxlastlickhit(i)=max(lastlickhit);
        end

        % plot heat map

        responsetrace10= nanmean(lickalign10);

        % plot average intensity trace
        xs10 = max(maxlastlickhit)-frameRate:(frameRate/secResolution):size(lickalign10,2);
        currxlim = [max(maxlastlickhit)-1.5*frameRate, max(maxlastlickhit)+4.5*frameRate];
        if loopcount == 1
            currxlim = [max(maxlastlickhit)-timeFrame, max(maxlastlickhit)+2*frameRate];
        end
        responsetrace10view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace10(abs(ceil(currxlim(1))):ceil(currxlim(2))));
        responsetrace10view(1:abs(floor(currxlim(1)))) = nan;
        
        hitAlignLL = responsetrace10view(abs(ceil(currxlim(1))):ceil(currxlim(2)));
 
        LLHitAlign = max(maxlastlickhit) * ones(1,2);
  

        %% FA trials aligned to last lick
        if FASize ~= 0
            rndlastlickFA=ceil(max(lastlickFA));
            lickalign11=NaN(length(b{1}.falseAlarmTrialNums(faDel:FASize)),100+rndlastlickFA);
            for i=1:length(b{1}.falseAlarmTrialNums(faDel:FASize)) % align licks
                tSize11 = length(hmuseFA((b{1}.falseAlarmTrialNums(i+faDel-1)-delTrial),:))-1;
                shift11 = ceil(max(lastlickFA) - lastlickFA(i)) + 1;
                lickalign11(i,shift11:(shift11+tSize11))=hmFALL(i,:);
            end
            for i=1:length(lastlickFA)
                maxlastlickFA(i)=max(lastlickFA);
            end

            % plot heat map

            responsetrace11= nanmean(lickalign11);

            % plot average intensity trace
            xs11 = max(maxlastlickFA)-frameRate:(frameRate/secResolution):size(lickalign11,2);    
            currxlim = [max(maxlastlickFA)-1.5*frameRate, max(maxlastlickFA)+4.5*frameRate];
            if loopcount == 1
                currxlim = [max(maxlastlickFA)-timeFrame, max(maxlastlickFA)+2*frameRate];
            end
            responsetrace11view = horzcat(zeros(1,abs(floor(currxlim(1)))),responsetrace11(abs(ceil(currxlim(1))):ceil(currxlim(2))));
            responsetrace11view(1:abs(floor(currxlim(1)))) = nan;
            
            FAAlignLL = responsetrace11(abs(ceil(currxlim(1))):ceil(currxlim(2)));

            LLFAAlign = max(maxlastlickFA) * ones(1,2);
        end

        %% Save data into cell
        if (loopcount == 1)
            timeNormData(sessionLoop,1) = mat2cell(hmhit,size(hmhit,1));
            timeNormData(sessionLoop,2) = mat2cell(hmhitLick,size(hmhitLick,1));
            timeNormData(sessionLoop,3) = mat2cell(hmhitpOnset,size(hmhitpOnset,1));
            if missSize ~= 0
                timeNormData(sessionLoop,4) = mat2cell(hmmiss,size(hmmiss,1));
                timeNormData(sessionLoop,5) = mat2cell(hmmisspOnset,size(hmmisspOnset,1));
            end
            if CRSize ~= 0                
                timeNormData(sessionLoop,6) = mat2cell(hmCR,size(hmCR,1));
                timeNormData(sessionLoop,7) = mat2cell(hmCRpOnset,size(hmCRpOnset,1));
                test = exist('CRAlignFL', 'var');
                if (test == 1)
                    timeNormData(sessionLoop,15) = mat2cell(CRAlignFL,size(CRAlignFL,1));
                end
                timeNormData(sessionLoop,20) = mat2cell(CRAlignPD,size(CRAlignPD,1));
            end
            if FASize ~= 0
                timeNormData(sessionLoop,8) = mat2cell(hmFA,size(hmFA,1));
                timeNormData(sessionLoop,9) = mat2cell(hmFApOnset,size(hmFApOnset,1));
                test = exist('FAAlignFAL', 'var');
                if (test == 1)
                    timeNormData(sessionLoop,11) = mat2cell(FAAlignFAL,size(FAAlignFAL,1));
                end
%                 timeNormData(sessionLoop,11) = mat2cell(FAAlignFAL,size(FAAlignFAL,1));
                if (test == 1)
                    timeNormData(sessionLoop,13) = mat2cell(FAAlignFL,size(FAAlignFL,1));
                end
%                 timeNormData(sessionLoop,13) = mat2cell(FAAlignFL,size(FAAlignFL,1));
                timeNormData(sessionLoop,18) = mat2cell(FAAlignPD,size(FAAlignPD,1));
                if (test == 1)
                    timeNormData(sessionLoop,22) = mat2cell(FAAlignLL,size(FAAlignLL,1));
                end                
%                 timeNormData(sessionLoop,22) = mat2cell(FAAlignLL,size(FAAlignLL,1));
            end
            
            timeNormData(sessionLoop,10) = mat2cell(hitAlignFAL,size(hitAlignFAL,1));
            timeNormData(sessionLoop,12) = mat2cell(hitAlignFL,size(hitAlignFL,1));
            test = exist('missAlignFL', 'var');
            if missSize ~= 0
                if (test == 1)
                    timeNormData(sessionLoop,14) = mat2cell(missAlignFL,size(missAlignFL,1));
                end
                timeNormData(sessionLoop,19) = mat2cell(missAlignPD,size(missAlignPD,1));
            end            
             timeNormData(sessionLoop,16) = mat2cell(hitAlignFLtrue,size(hitAlignFLtrue,1));
            timeNormData(sessionLoop,17) = mat2cell(hitAlignPD,size(hitAlignPD,1));
            timeNormData(sessionLoop,21) = mat2cell(hitAlignLL,size(hitAlignLL,1));
            timeNormData(sessionLoop,23) = mat2cell(onset,size(onset,1));

        elseif loopcount == 2
            frameNormData(sessionLoop,1) = mat2cell(hmhit,size(hmhit,1));
            frameNormData(sessionLoop,2) = mat2cell(hmhitLick,size(hmhitLick,1));
            frameNormData(sessionLoop,3) = mat2cell(hmhitpOnset,size(hmhitpOnset,1));
            if missSize ~= 0
                frameNormData(sessionLoop,4) = mat2cell(hmmiss,size(hmmiss,1));
                frameNormData(sessionLoop,5) = mat2cell(hmmisspOnset,size(hmmisspOnset,1));
            end
            if FASize ~= 0
                frameNormData(sessionLoop,8) = mat2cell(hmFA,size(hmFA,1));
                frameNormData(sessionLoop,9) = mat2cell(hmFApOnset,size(hmFApOnset,1));
                test = exist('FAAlignFAL', 'var');
                if (test == 1)
                    frameNormData(sessionLoop,11) = mat2cell(FAAlignFAL,size(FAAlignFAL,1));
                end
%                 frameNormData(sessionLoop,11) = mat2cell(FAAlignFAL,size(FAAlignFAL,1));
                test = exist('FAAlignFL', 'var');
                if (test == 1)
                     frameNormData(sessionLoop,13) = mat2cell(FAAlignFL,size(FAAlignFL,1));
                end
%                 frameNormData(sessionLoop,13) = mat2cell(FAAlignFL,size(FAAlignFL,1));
                frameNormData(sessionLoop,18) = mat2cell(FAAlignPD,size(FAAlignPD,1));
                test = exist('FAAlignLL', 'var');
                if (test == 1)
                     frameNormData(sessionLoop,22) = mat2cell(FAAlignLL,size(FAAlignLL,1));
                end                
%                 frameNormData(sessionLoop,22) = mat2cell(FAAlignLL,size(FAAlignLL,1));
            end
            if CRSize ~= 0
                frameNormData(sessionLoop,6) = mat2cell(hmCR,size(hmCR,1));
                frameNormData(sessionLoop,7) = mat2cell(hmCRpOnset,size(hmCRpOnset,1));
                test = exist('CRAlignFL', 'var');
                if (test == 1)
                    frameNormData(sessionLoop,15) = mat2cell(CRAlignFL,size(CRAlignFL,1));
                end
                frameNormData(sessionLoop,20) = mat2cell(CRAlignPD,size(CRAlignPD,1));
            end
            frameNormData(sessionLoop,10) = mat2cell(hitAlignFAL,size(hitAlignFAL,1));
            frameNormData(sessionLoop,12) = mat2cell(hitAlignFL,size(hitAlignFL,1));
            test = exist('missAlignFL', 'var');
            if missSize ~= 0
                if (test == 1)
                    frameNormData(sessionLoop,14) = mat2cell(missAlignFL,size(missAlignFL,1));
                end
                frameNormData(sessionLoop,19) = mat2cell(missAlignPD,size(missAlignPD,1));
            end            
            frameNormData(sessionLoop,16) = mat2cell(hitAlignFLtrue,size(hitAlignFLtrue,1));
            frameNormData(sessionLoop,17) = mat2cell(hitAlignPD,size(hitAlignPD,1));
            frameNormData(sessionLoop,21) = mat2cell(hitAlignLL,size(hitAlignLL,1));
            frameNormData(sessionLoop,23) = mat2cell(onset,size(onset,1));
        elseif loopcount == 3
            baseData(sessionLoop,1) = mat2cell(hmhit,size(hmhit,1));
            baseData(sessionLoop,2) = mat2cell(hmhitLick,size(hmhitLick,1));
            baseData(sessionLoop,3) = mat2cell(hmhitpOnset,size(hmhitpOnset,1));
            if missSize ~= 0
                baseData(sessionLoop,4) = mat2cell(hmmiss,size(hmmiss,1));
                baseData(sessionLoop,5) = mat2cell(hmmisspOnset,size(hmmisspOnset,1));
            end
            if FASize ~= 0
                baseData(sessionLoop,8) = mat2cell(hmFA,size(hmFA,1));
                baseData(sessionLoop,9) = mat2cell(hmFApOnset,size(hmFApOnset,1));
                test = exist('FAAlignFAL', 'var');
                if (test == 1)
                     baseData(sessionLoop,11) = mat2cell(FAAlignFAL,size(FAAlignFAL,1));
                end                
%                 baseData(sessionLoop,11) = mat2cell(FAAlignFAL,size(FAAlignFAL,1));
                test = exist('FAAlignFL', 'var');
                if (test == 1)
                     baseData(sessionLoop,13) = mat2cell(FAAlignFL,size(FAAlignFL,1));
                end     
%                 baseData(sessionLoop,13) = mat2cell(FAAlignFL,size(FAAlignFL,1));
                baseData(sessionLoop,18) = mat2cell(FAAlignPD,size(FAAlignPD,1));
                test = exist('FAAlignLL', 'var');
                if (test == 1)
                     baseData(sessionLoop,22) = mat2cell(FAAlignLL,size(FAAlignLL,1));
                end                   
%                 baseData(sessionLoop,22) = mat2cell(FAAlignLL,size(FAAlignLL,1));
            end
            if CRSize ~= 0
                baseData(sessionLoop,6) = mat2cell(hmCR,size(hmCR,1));
                baseData(sessionLoop,7) = mat2cell(hmCRpOnset,size(hmCRpOnset,1));
                test = exist('CRAlignFL', 'var');
                if (test == 1)
                    baseData(sessionLoop,15) = mat2cell(CRAlignFL,size(CRAlignFL,1));
                end
                baseData(sessionLoop,20) = mat2cell(CRAlignPD,size(CRAlignPD,1));
            end
            baseData(sessionLoop,10) = mat2cell(hitAlignFAL,size(hitAlignFAL,1));
            baseData(sessionLoop,12) = mat2cell(hitAlignFL,size(hitAlignFL,1));
            test = exist('missAlignFL', 'var');
            if missSize ~= 0
                if (test == 1)
                    baseData(sessionLoop,14) = mat2cell(missAlignFL,size(missAlignFL,1));
                end
                baseData(sessionLoop,19) = mat2cell(missAlignPD,size(missAlignPD,1));
            end
            baseData(sessionLoop,16) = mat2cell(hitAlignFLtrue,size(hitAlignFLtrue,1));
            baseData(sessionLoop,17) = mat2cell(hitAlignPD,size(hitAlignPD,1));
            baseData(sessionLoop,21) = mat2cell(hitAlignLL,size(hitAlignLL,1));
            baseData(sessionLoop,23) = mat2cell(onset,size(onset,1));
        end
    end
    clearvars -except fileNames baseData timeNormData frameNormData sessionLoop
    cd('..');
    close all
    toc
end

close all
comboSave = "";
for i = 1:length(fileNames)
    fname = char(fileNames(i));
    comboSave = strcat(comboSave,fname(6:9),'-');
end
comboSave = char(comboSave);
comboSave = comboSave(1:end-1);
comboSave = strcat(comboSave,'.mat');
save(comboSave);
disp('dbcont to continue, dbquit to stop');
% keyboard
% Color_Plots;