%Krista Wahlstrom, 2023
% This is the behavioral analysis script for BLAES Aim 1 Long Delay for the Test
%Phase (retreival) session. It includes a correction for Hit Rates and
%False Alarm rates that are 0 or 1 using the Loglinear (Hautus 1995)
%approach (0.5 added to the number of hits and 0.5 added to the number of false alarms 
% and 1 added to both the # of signal trials and # of noise trials. These corrections are
%made REGARDLESS if extreme 0 or 1 values are present for the False Alarm
%and Hit Rates). The corrections are made within each individual Figure (for the dprime figures only. all other figures reflect true values),
%and False Alarm occurrences/Hit occurrences are not altered in the Figure
%legends. The legends still reflect the actual # of FAs or Hits the
%participant had.

function BLAES_Aim1_TestPhaseBehavioralAnalysis_Loglinear_Correction()

%% StimulusCodes
% 1 - 1000 Images
% 1001 - 2000 Image onscreen for 1s between Yes/No response and Sure/Not Sure response (only appplies to the Test phase)
% 2501-2504 Instructions screens
% 2601 Sync pulse
% 2701 No Stim (only applies to the STUDY PHASE)
% 2702 Stim (only applies to the STUDY PHASE)
% 5200 ISI (w/fixation cross onscreen)
% 2801 - 5115 Variable timed fixation cross
% 5301 - 6319 Confidence ratings (only applies to the TEST PHASE)


%% Load data
clear;
close all;

addpath(genpath(fullfile(cd,'BCI2000Tools')))

subjID         = 'BJH029';
RetrievalDelay = 'One-Day'; %Set to One-Day or Long-Delay depending on which retrieval test you want to analyze

%Edit these to include all 4 response options or only the 2 "sure" options
NoResponse     = [67];
YesResponse    = [78];
ResponseKeys   = [NoResponse, YesResponse]; % 67 sure no, 86 not-sure no, 66 not-sure yes, 78 sure yes

if size(ResponseKeys,2) == 2
    SureResponseString = '_Sure';
else
    SureResponseString = '';
end



LoadStoredData = 1;
figsize        = [100 100 1200 800];


      
%% Get relevant .dat files
d = dir(fullfile(cd,'data',subjID,'Test', RetrievalDelay,'*.dat'));



%% Extract behavioral data
if LoadStoredData && exist(fullfile(cd,'data',subjID,'Test',RetrievalDelay,strcat(subjID,'_Test_Data_',RetrievalDelay,SureResponseString,'_Loglinear','.mat')))
    load(fullfile(cd,'data',subjID,'Test',RetrievalDelay,strcat(subjID,'_Test_Data_',RetrievalDelay,SureResponseString,'_Loglinear','.mat')))
else
    iter2 = 1;
    for file = 1:size(d,1)
        
        [~, states, parameters] = load_bcidat(fullfile(d(file).folder,d(file).name));
        pause(1);
        
        seq         = parameters.Sequence.NumericValue;
        seqResponse = seq;
        KD          = states.KeyDown;
        StimCode    = states.StimulusCode;


        
        %% clean up sequence
        % select only image stimuli
        seq(seq<1) = [];
        seq(seq>1000) = [];
        
        %confidence ratings stimulus codes
        seqResponse(seqResponse<5301)  = [];
        seqResponse(seqResponse>6319) = [];


        %Create a copy of the keypresses file
        KDCopy = KD;
        KDCopy(KDCopy==0) = [];


        %Create keyPresses variable which lists the keypresses, the stimcode, and the stimcode index at which that key press occurred 
        cnt = 1;
        keyPresses = zeros(length(KDCopy), 3);
        for i = 1:length(KD)
            if(KD(i) ~= 0)
                keyPresses(cnt,1) = KD(i);
                keyPresses(cnt,2) = StimCode(i);
                keyPresses(cnt,3) = i;
                cnt = cnt + 1;
            end
        end
        
        
        
     


        
        %% Compile data into single matrix
        for i = 1:length(seq)
            collectData{iter2,1} = parameters.Stimuli.Value{2,seq(i)};          % picture filename
            collectData{iter2,2} = str2num(parameters.Stimuli.Value{8,seq(i)}); % Stim or no stim during encoding
            collectData{iter2,3} = parameters.Stimuli.Value{9,seq(i)};          % target (old) or new
            collectData{iter2,4} = seq(i);                                      % stimulus code for image
            
            %Collected 'Yes'/'No' key presses from when the image stimulus
            %code was present
            idx = ismember(keyPresses(:,2), seq(i));
            keyPressesForSeq = keyPresses(idx, :);
            
            %Collect the first valid key press (if patient presses multiple
            %keys during the image, take the first valid (left or right
            %arrow) press)
            pressForImage = 0;
            for j = 1: size(keyPressesForSeq)
                if(pressForImage == 0)
                    pressForImage = keyPressesForSeq(j, 1);
                elseif(keyPressesForSeq(j, 1) == 37 || keyPressesForSeq(j, 1) == 39)
                    pressForImage = keyPressesForSeq(j, 1);
                end
    
                if(pressForImage == 37 || pressForImage == 39)
                    break;
                end
            end

            if pressForImage == 37
                collectData{iter2,5} = 'Yes';
            elseif pressForImage == 39
                collectData{iter2,5} = 'No';
            else
                collectData{iter2,5} = 'Non-Response Key';
            end

            collectData{iter2, 6} = pressForImage;
            



            %Collected 'Sure'/'Not Sure' key presses from when the confidence rating stimulus
            %code was present
            idx_confidence = ismember(keyPresses(:,2), seqResponse(i));
            keyPressesForConf = keyPresses(idx_confidence, :);
            
            %Collect the first valid key press (if patient presses multiple
            %keys during the image/confidence rating period, take the first valid (left or right
            %arrow) press)
            pressForConf = 0;
            for j = 1: size(keyPressesForConf)
                if(pressForConf == 0)
                    pressForConf = keyPressesForConf(j, 1);
                elseif(keyPressesForConf(j, 1) == 37 || keyPressesForConf(j, 1) == 39)
                    pressForConf = keyPressesForConf(j, 1);
                end
    
                if(pressForConf == 37 || pressForConf == 39)
                    break;
                end
            end

            if pressForConf == 37
                collectData{iter2,7} = 'Sure';
            elseif pressForConf == 39
                collectData{iter2,7} = 'Not-Sure';
            else
                collectData{iter2,7} = 'Non-Response Key';
            end

            collectData{iter2, 8} = pressForConf;
            

            %Code the Yes/No and Sure/NotSure responses the same as in
            %BLAES Aim2.1, where 67 = sure no, 86 = not-sure no, 66 =
            %not-sure yes, 78 = sure yes (so we have one value that
            %encompasses both key presses in this two-stage paradigm)
            if (collectData {iter2,6} == 39 && collectData {iter2,8} == 37)
                collectData {iter2,9} = 67;
            elseif (collectData {iter2,6} == 39 && collectData {iter2,8} == 39)
                collectData {iter2,9} = 86;
            elseif (collectData {iter2,6} == 37 && collectData {iter2,8} == 39)
                collectData {iter2,9} = 66;
            elseif (collectData {iter2,6} == 37 && collectData {iter2,8} == 37)
                collectData {iter2,9} = 78;
            end
            




            iter2 = iter2 + 1;
        end
        
    end
    
    save(fullfile(cd,'data',subjID,'Test',RetrievalDelay,strcat(subjID,'_Test_Data_',RetrievalDelay,SureResponseString,'_Loglinear','.mat')),'collectData')
end


%% Behavioral Analysis
% hit/miss/false alarm/correct rejection
% Yes (targ/old) No (New)       Sure           Not Sure
% KeyDown == 37  KeyDown == 39  KeyDown == 37  KeyDown == 39

%Response values that encompass both Yes/No and Sure/NotSure
% 67 = sure no, 86 = not-sure no, 66 = not-sure yes, 78 = sure yes



% Get total number of possible HITS for each category (target/old, object/image, stim/nostim) ---total possible
% HITS and total possible MISSES are the same;
TotalObjectHitStim    = GetBA_MaxPossible(collectData, 'Targ',  1);
TotalObjectHitNoStim  = GetBA_MaxPossible(collectData, 'Targ',  0);


% Get total number of possible FALSE ALARMS for each category
% (new no stim = 0 (new images are always no stim)
% Total number of possible false alarms and rejections is the same
TotalObjectFA  = GetBA_MaxPossible(collectData, 'New', 0);



for i = 1:size(collectData,1)
    
    % Get counts for all patient responses (hits, misses, false alarms,
    % rejections)
    
    ObjectHit_Stim(i,1)     = CheckBAConditions(collectData(i,:), 'Targ',  1, YesResponse);
    ObjectHit_NoStim(i,1)   = CheckBAConditions(collectData(i,:), 'Targ',  0, YesResponse);
    
    
    ObjectMiss_Stim(i,1)    = CheckBAConditions(collectData(i,:), 'Targ',  1, NoResponse);
    ObjectMiss_NoStim(i,1)  = CheckBAConditions(collectData(i,:), 'Targ',  0, NoResponse);
    
    
    ObjectFalseAlarm(i,1)   = CheckBAConditions(collectData(i,:), 'New',  0, YesResponse);
    
    ObjectRejection(i,1)    = CheckBAConditions(collectData(i,:), 'New',  0, NoResponse);
    
end


%% Plot Figures


%%  Total hits across all categories
% combining stim and nostim hits for object images
labelstring = {'NoStim + Stim'};
textStartY = 0.75;
textStepY  = 0.05;

TrainedImagesCombinedStimFig = figure('Position',figsize);
b = bar((sum(ObjectHit_NoStim + ObjectHit_Stim))/...
    (TotalObjectHitNoStim + TotalObjectHitStim));
b(1).FaceColor = [0.4940 0.1840 0.5560]; % purple
text(1.55,textStartY,'number of hits:','FontSize',18,'fontweight','bold')
text(1.55,textStartY-1*textStepY,labelstring,'FontSize',18)
text(1.8,textStartY-1*textStepY,num2str(sum(ObjectHit_NoStim + ObjectHit_Stim)),'FontSize',18)
ylabel('Occurrences (% of Total)')
set(gca,'YTick',0:0.25:1,'YTickLabel',100*round(0:0.25:1,2))
set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
title([subjID, ' ', RetrievalDelay ' Hit Rate - Combined Stim/NoStim'])
set(gca,'XTick',1,'XTickLabel',{'Stim/NoStim'})
axis([0.5 2.0 0 1])

fprintf('Saving trained images - combined stim figure...\n')
savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_TrainedImagesCombinedStim.png'));

saveas(TrainedImagesCombinedStimFig, [savefile(1:end-4) '.png'],'png');
fprintf('Done\n')

%% Hit rate for objects no stim, objects stim
% break out stim and nostim (hit rate out of total images in each category)

labelstring = {'NoStim','Stim'};
textStartY = 0.75;
textStepY  = 0.05;

TrainedImagesFig = figure('Position',figsize);
b = bar([((sum(ObjectHit_NoStim))/(TotalObjectHitNoStim));    ((sum(ObjectHit_Stim))/(TotalObjectHitStim))],'FaceColor','Flat');

b.CData(1,:) = [0 0 1]; %blue
b.CData(2,:) = [1 0 0]; %red



%Create hatches for bars
drawnow
vertdat = b.Face.VertexData(:,1:4);

hp = patch(vertdat(1,:),vertdat(2,:),[1 1 1]);

hhf = hatchfill2(hp,'single','HatchAngle',45,'HatchDensity',60,'hatchcolor',[0 0 1]);

vertdat = b.Face.VertexData(:,5:8);

hp = patch(vertdat(1,:),vertdat(2,:),[1 1 1]);

hhf = hatchfill2(hp,'cross','HatchAngle',45,'HatchDensity',40,'hatchcolor',[1 0 0]);




text(2.45,textStartY,'number of hits:','FontSize',18,'fontweight','bold')
for i = 1:size(labelstring,2)
    text(2.45,textStartY-i*textStepY,labelstring(i),'FontSize',18)
end
text(2.8,textStartY-1*textStepY,num2str(sum(ObjectHit_NoStim)),'FontSize',18)
text(2.8,textStartY-2*textStepY,num2str(sum(ObjectHit_Stim)),'FontSize',18)

set(gca,'XTick',[1 2],'XTickLabel',{'NoStim','Stim'})
ylabel('Occurrences (% of Total)')
set(gca,'YTick',0:0.25:1,'YTickLabel',100*round(0:0.25:1,2))
set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
axis([0.5 3.0 0 1])
title([subjID, ' ', RetrievalDelay ' Hit Rate (out of total images)'])

% %Create Legend
% legendData = {'No Stim', 'Stim'};
% [legend_h, object_h, plot_h, text_str] = legendflex(b, legendData, 'Padding', [2, 2, 10], 'FontSize', 18, 'Location', 'NorthEast');
% % object_h(1) is the first bar's text
% % object_h(2) is the second bar's text
% % object_h(3) is the first bar's patch
% % object_h(4) is the second bar's patch
% %
% % Set the two patches within the legend
% hatchfill2(object_h(3), 'single','HatchAngle',45, 'HatchDensity',60/4,'hatchcolor',[0 0 1]);
% hatchfill2(object_h(4), 'cross','HatchAngle',45,'HatchDensity',40/4,'hatchcolor',[1 0 0]);


fprintf('Saving trained images figure...\n')
savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_TrainedImages.png'));

saveas(TrainedImagesFig, [savefile(1:end-4) '.png'],'png');
fprintf('Done\n')


%% Hit rate for objects no stim, objects stim 
% (ACTUAL hit rate -- Sure only for total number of "sure" images)

%This analysis/figure only runs when the two "sure" key responses are given
%in lines 36/37
if size(ResponseKeys,2) == 2

    labelstring = {'NoStim','Stim'};
    textStartY = 0.75;
    textStepY  = 0.05;
    
    ActualHitRateSureFig = figure('Position',figsize);
    b = bar([((sum(ObjectHit_NoStim))/(sum(ObjectHit_NoStim) + sum(ObjectMiss_NoStim)));    ((sum(ObjectHit_Stim))/(sum(ObjectHit_Stim) + sum(ObjectMiss_Stim)))],'FaceColor','Flat');

    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];

    
    %Create hatches for bars
    drawnow
    vertdat = b.Face.VertexData(:,1:4);

    hp = patch(vertdat(1,:),vertdat(2,:),[1 1 1]);

    hhf = hatchfill2(hp,'single','HatchAngle',45,'HatchDensity',60,'hatchcolor',[0 0 1]);

    vertdat = b.Face.VertexData(:,5:8);

    hp = patch(vertdat(1,:),vertdat(2,:),[1 1 1]);

    hhf = hatchfill2(hp,'cross','HatchAngle',45,'HatchDensity',40,'hatchcolor',[1 0 0]);



    
    text(2.45,textStartY,'number of hits:','FontSize',18,'fontweight','bold')
    for i = 1:size(labelstring,2)
        text(2.45,textStartY-i*textStepY,labelstring(i),'FontSize',18)
    end
  
    
    text(2.8,textStartY-1*textStepY,strcat(num2str(sum(ObjectHit_NoStim)),' (out of ',num2str(sum(ObjectHit_NoStim) + sum(ObjectMiss_NoStim)), ')'),'FontSize',18)
    text(2.8,textStartY-2*textStepY,strcat(num2str(sum(ObjectHit_Stim)),' (out of ',num2str(sum(ObjectHit_Stim) + sum(ObjectMiss_Stim)), ')'),'FontSize',18)

    set(gca,'XTick',[1 2],'XTickLabel',{'NoStim','Stim'})
    ylabel('Occurrences (% of Sure)')
    set(gca,'YTick',0:0.25:1,'YTickLabel',100*round(0:0.25:1,2))
    set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
    axis([0.5 3.0 0 1])
    title([subjID, ' ', RetrievalDelay ' Actual Hit Rate (out of sure only responses)'])
    
%     %Create Legend
%     legendData = {'No Stim', 'Stim'};
%     [legend_h, object_h, plot_h, text_str] = legendflex(b, legendData, 'Padding', [2, 2, 10], 'FontSize', 18, 'Location', 'NorthEast');
%     % object_h(1) is the first bar's text
%     % object_h(2) is the second bar's text
%     % object_h(3) is the first bar's patch
%     % object_h(4) is the second bar's patch
%     %
%     % Set the two patches within the legend
%     hatchfill2(object_h(3), 'single','HatchAngle',45, 'HatchDensity',60/4,'hatchcolor',[0 0 1]);
%     hatchfill2(object_h(4), 'cross','HatchAngle',45,'HatchDensity',40/4,'hatchcolor',[1 0 0]);
    
    
    fprintf('Saving trained images figure...\n')
    savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_ActualHitRateSURE.png'));
    
    saveas(ActualHitRateSureFig, [savefile(1:end-4) '.png'],'png');
    fprintf('Done\n')

end


%% False Alarm Rate
%(out of total images shown for each category)
labelstring = {'Objects'};
textStartY = 0.75;
textStepY  = 0.05;

NovelImagesFig = figure('Position',figsize);
b = bar(((sum(ObjectFalseAlarm))/(TotalObjectFA)),'FaceColor','Flat');
b(1).FaceColor = [0.5 0.5 0.5];
text(1.45,textStartY,'number of FAs:','FontSize',18,'fontweight','bold')
for i = 1:size(labelstring,2)
    text(1.45,textStartY-i*textStepY,labelstring(i),'FontSize',18)
end
text(1.85,textStartY-1*textStepY,num2str(sum(ObjectFalseAlarm)),'FontSize',18)
set(gca,'XTick',1,'XTickLabel',{'Objects'})
ylabel('Occurrences (% of Total)')
set(gca,'YTick',0:0.25:1,'YTickLabel',100*round(0:0.25:1,2))
set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
axis([0.5 3.0 0 1])
title([subjID, ' ', RetrievalDelay ' False Alarm (out of total images)'])

fprintf('Saving novel images Figure...\n')
savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_NovelImages.png'));

saveas(NovelImagesFig, [savefile(1:end-4) '.png'],'png');
fprintf('Done\n')


%% False alarm rate (ACTUAL FA rate -- Sure only for total number of "sure" images)

%This analysis/figure only runs when the two "sure" key responses are given
%in lines 36/37
if size(ResponseKeys,2) == 2

    labelstring = {'Objects'};
    textStartY = 0.75;
    textStepY  = 0.05;
    
    ActualFARateSureFig = figure('Position',figsize);
    b = bar(((sum(ObjectFalseAlarm))/(sum(ObjectFalseAlarm) + sum(ObjectRejection))),'FaceColor','Flat');
    b(1).FaceColor = [0.5 0.5 0.5];
    text(1.45,textStartY,'number of FAs:','FontSize',18,'fontweight','bold')
    for i = 1:size(labelstring,2)
        text(1.45,textStartY-i*textStepY,labelstring(i),'FontSize',18)
    end
    text(1.85,textStartY-1*textStepY,strcat(num2str(sum(ObjectFalseAlarm)),' (out of ',num2str(sum(ObjectFalseAlarm) + sum(ObjectRejection)), ')'),'FontSize',18)

    set(gca,'XTick',1,'XTickLabel',{'Objects'})
    ylabel('Occurrences (% of Sure)')
    set(gca,'YTick',0:0.25:1,'YTickLabel',100*round(0:0.25:1,2))
    set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
    axis([0.5 3.0 0 1])
    title([subjID, ' ', RetrievalDelay ' Actual False Alarm (out of sure only responses)'])
    
    fprintf('Saving novel images Figure...\n')
    savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_ActualFArateSure.png'));
    
    saveas(ActualFARateSureFig, [savefile(1:end-4) '.png'],'png');
    fprintf('Done\n')

end

%% d prime 
% Breakout stim, no stim (includes loglinear correction)
[dp(1,1), c(1)] = dprime_simple((sum(ObjectHit_NoStim) + 0.5)/(TotalObjectHitNoStim + 1),(sum(ObjectFalseAlarm) + 0.5)/(TotalObjectFA + 1));         % nostim 
[dp(1,2), c(2)] = dprime_simple((sum(ObjectHit_Stim) + 0.5)/(TotalObjectHitStim + 1),(sum(ObjectFalseAlarm) + 0.5)/(TotalObjectFA + 1));             % stim 

dp(isinf(dp)) = NaN;

%For y axis limit locked at 4.0
ylim = [min([min(dp,[],'omitnan'),0])-0.1 max(4.0)];

%For variable y axis limit
%ylim = [min([min(dp,[],'omitnan'),0])-0.1 max([max(dp,[],'omitnan'),0])+0.1];

          
clabelstring    = {'objects nostim', 'objects stim'};
textStartY = max(max(dp,[],'omitnan')) - 0.15*max(max(dp,[],'omitnan'));
textStepY  = max(max(dp,[],'omitnan'))/20;

dprimefig = figure('Position',figsize);
b = bar(dp,'FaceColor','flat');

b.CData(1,:) = [0 0 1];
b.CData(2,:) = [1 0 0];



text(3.0,3.0,'Criterion:','FontSize',18,'fontweight','bold')
text(3.0,2.3,'Dprime:','FontSize',18,'fontweight','bold')

%Plot dprime values
DText1 = {'NoStim: '};
DText1a = {strcat('  ',num2str(dp(1,1),3))};
DText2 = {'Stim: '};
DText2a = {strcat('  ',num2str(dp(1,2),3))};


text(3.0, 2.2, DText1,'FontSize',15)
text(3.5, 2.2, DText1a,'FontSize',15)
text(3.0, 2.1, DText2,'FontSize',15)
text(3.5, 2.1, DText2a,'FontSize',15)






%Plot Criterion text for when yaxis is locked at dprime of 4.0
Text1 = {'NoStim: '};
Text1a = {strcat('  ',num2str(c(1),3))};
Text2 = {'Stim: '};
Text2a = {strcat('  ',num2str(c(2),3))};



text(3.0, 2.9, Text1,'FontSize',15)
text(3.5, 2.9, Text1a,'FontSize',15)
text(3.0, 2.8, Text2,'FontSize',15)
text(3.5, 2.8, Text2a,'FontSize',15)


%Plot criterion text if you have yaxis changing based on dprime vaues and
%not locked
% for i = 1:numel(c)
%     text(3.35,textStartY-i*textStepY,clabelstring(i),'FontSize',18)
%     if c(i)<0
%     	text(3.965,textStartY-i*textStepY,num2str(c(i),3),'FontSize',18)
%     else
%         text(4,textStartY-i*textStepY,num2str(c(i),3),'FontSize',18)
%     end
% end


set(gca,'XTick',[1 2],'XTickLabel',{'NoStim','Stim'})
ylabel('dprime')
set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
axis([0.5 4.5 ylim])
title([subjID, ' ', RetrievalDelay ' Discrimination Index (out of total images)'])

%Create plot text that lists the difference in dprime for stim vs. nostim
plottext1 = {'Dprime Diff: '};
plottext2 = {strcat('  ',num2str(abs(dprime_simple((sum(ObjectHit_NoStim) + 0.5)/(TotalObjectHitNoStim + 1),(sum(ObjectFalseAlarm) + 0.5)/(TotalObjectFA + 1))-dprime_simple((sum(ObjectHit_Stim) + 0.5)/(TotalObjectHitStim + 1),(sum(ObjectFalseAlarm) + 0.5)/(TotalObjectFA + 1))),3))};

%Set the positioning of the above text
text(0.8, 3.8, plottext1,'FontSize',19)
text(0.8, 3.6, plottext2,'FontSize',19)





fprintf('Saving d prime figure...\n')
savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_dprime.png'));

saveas(dprimefig, [savefile(1:end-4) '.png'],'png');
fprintf('Done\n')


%% d prime (ACTUAL dprime -- Sure only for total number of "sure" images)

%This analysis/figure only runs when the two "sure" key responses are given
%in lines 36/37
if size(ResponseKeys,2) == 2

    [dp(1,1), c(1)] = dprime_simple((sum(ObjectHit_NoStim) + 0.5)/(sum(ObjectHit_NoStim) + sum(ObjectMiss_NoStim) + 1),(sum(ObjectFalseAlarm) + 0.5) /(sum(ObjectFalseAlarm) + sum(ObjectRejection) + 1));         % nostim 
    [dp(1,2), c(2)] = dprime_simple((sum(ObjectHit_Stim) + 0.5)/(sum(ObjectHit_Stim) + sum(ObjectMiss_Stim) + 1),(sum(ObjectFalseAlarm) + 0.5)/(sum(ObjectFalseAlarm) + sum(ObjectRejection) + 1));             % stim 
    
    dp(isinf(dp)) = NaN;
    
    %For y axis limit locked at 4.0
    ylim = [min([min(dp,[],'omitnan'),0])-0.1 max(4.0)];
    
    %For variable y axis limit
    %ylim = [min([min(dp,[],'omitnan'),0])-0.1 max([max(dp,[],'omitnan'),0])+0.1];
    
              
    clabelstring    = {'objects nostim', 'objects stim'};
    textStartY = max(max(dp,[],'omitnan')) - 0.15*max(max(dp,[],'omitnan'));
    textStepY  = max(max(dp,[],'omitnan'))/20;
    
  
    Actualdprimefig = figure('Position',figsize);
    
    b = bar(dp,'FaceColor','flat');

    b.CData(1,:) = [0 0 1];
    b.CData(2,:) = [1 0 0];
    

    text(3.0,3.0,'Criterion:','FontSize',18,'fontweight','bold')
    text(3.0,2.3,'Dprime:','FontSize',18,'fontweight','bold')
    
    %Plot dprime values
    DText1 = {'NoStim: '};
    DText1a = {strcat('  ',num2str(dp(1,1),3))};
    DText2 = {'Stim: '};
    DText2a = {strcat('  ',num2str(dp(1,2),3))};
    
    
    text(3.0, 2.2, DText1,'FontSize',15)
    text(3.5, 2.2, DText1a,'FontSize',15)
    text(3.0, 2.1, DText2,'FontSize',15)
    text(3.5, 2.1, DText2a,'FontSize',15)
    
    
    
    
    
    
    %Plot Criterion text for when yaxis is locked at dprime of 4.0
    Text1 = {'NoStim: '};
    Text1a = {strcat('  ',num2str(c(1),3))};
    Text2 = {'Stim: '};
    Text2a = {strcat('  ',num2str(c(2),3))};
    
    
    
    text(3.0, 2.9, Text1,'FontSize',15)
    text(3.5, 2.9, Text1a,'FontSize',15)
    text(3.0, 2.8, Text2,'FontSize',15)
    text(3.5, 2.8, Text2a,'FontSize',15)
    
    
    %Plot criterion text if you have yaxis changing based on dprime vaues and
    %not locked
    % for i = 1:numel(c)
    %     text(3.35,textStartY-i*textStepY,clabelstring(i),'FontSize',18)
    %     if c(i)<0
    %     	text(3.965,textStartY-i*textStepY,num2str(c(i),3),'FontSize',18)
    %     else
    %         text(4,textStartY-i*textStepY,num2str(c(i),3),'FontSize',18)
    %     end
    % end
    
    
    set(gca,'XTick',[1 2],'XTickLabel',{'NoStim','Stim'})
    ylabel('dprime')
    set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
    axis([0.5 4.5 ylim])
    title([subjID, ' ', RetrievalDelay ' Actual Discrimination Index (out of sure only responses)'])
    
    %Create plot text that lists the difference in dprime for stim vs. nostim
    plottext1 = {'Dprime Diff: '};
    plottext2 = {strcat('  ',num2str(abs(dprime_simple((sum(ObjectHit_NoStim) + 0.5)/(sum(ObjectHit_NoStim) + sum(ObjectMiss_NoStim) + 1),(sum(ObjectFalseAlarm) + 0.5) /(sum(ObjectFalseAlarm) + sum(ObjectRejection) + 1))-dprime_simple((sum(ObjectHit_Stim) + 0.5)/(sum(ObjectHit_Stim) + sum(ObjectMiss_Stim) + 1),(sum(ObjectFalseAlarm) + 0.5)/(sum(ObjectFalseAlarm) + sum(ObjectRejection) + 1))),3))};
    
    %Set the positioning of the above text
    text(0.8, 3.8, plottext1,'FontSize',19)
    text(0.8, 3.6, plottext2,'FontSize',19)
    
    
    
    
    
    fprintf('Saving d prime figure...\n')
    savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_ActualSuredprime.png'));
    
    saveas(Actualdprimefig, [savefile(1:end-4) '.png'],'png');
    fprintf('Done\n')

end



%% d prime - combined stim and no stim
[dp_CombinedStim(1,1), c_CombinedStim(1)] = dprime_simple((sum(ObjectHit_NoStim + ObjectHit_Stim) + 0.5)/...
                                                (TotalObjectHitNoStim + TotalObjectHitStim + 1),(sum(ObjectFalseAlarm) + 0.5)/(TotalObjectFA + 1));           % item
                                  
dp_CombinedStim(isinf(dp_CombinedStim)) = NaN;

ylim = [min([min(dp_CombinedStim,[],'omitnan'),0])-0.1 max([max(dp_CombinedStim,[],'omitnan'),0])+0.1];
                                                
clabelstring    = {'Objects'};
textStartY = max(max(dp_CombinedStim,[],'omitnan')) - 0.15*max(max(dp_CombinedStim,[],'omitnan'));
textStepY  = max(max(dp_CombinedStim,[],'omitnan'))/20;

dprimeCombinedStimfig = figure('Position',figsize);
b = bar(dp_CombinedStim);
b(1).FaceColor = [0.3010 0.7450 0.9330]; % blue
text(2.0,textStartY,'Criterion:','FontSize',18,'fontweight','bold')
for i = 1:numel(c_CombinedStim)
    text(2.0,textStartY-i*textStepY,clabelstring(i),'FontSize',18)
    if c_CombinedStim(i)<0
    	text(3.965,textStartY-i*textStepY,num2str(c_CombinedStim(i)),'FontSize',18)
    else
        text(4,textStartY-i*textStepY,num2str(c_CombinedStim(i)),'FontSize',18)
    end
end
set(gca,'XTick',1,'XTickLabel',{'Objects'})
ylabel('dprime')
set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
axis([0.5 4.5 ylim])
title([subjID, ' ', RetrievalDelay ' Discrimination Index - Combined Stim/NoStim (out of total images)'])

fprintf('Saving d prime with combined stim figure...\n')
savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_dprimeCombinedStim.png'));

saveas(dprimeCombinedStimfig, [savefile(1:end-4) '.png'],'png');
fprintf('Done\n')



%% d prime - combined stim and no stim (ACTUAL dprime -- Sure only for total number of "sure" images)

%This analysis/figure only runs when the two "sure" key responses are given
%in lines 36/37
if size(ResponseKeys,2) == 2

    [dp_CombinedStim(1,1), c_CombinedStim(1)] = dprime_simple((sum(ObjectHit_NoStim + ObjectHit_Stim) + 0.5)/...
                                                    (sum(ObjectHit_NoStim + ObjectHit_Stim) + sum(ObjectMiss_NoStim + ObjectMiss_Stim) + 1),(sum(ObjectFalseAlarm) + 0.5)/(sum(ObjectFalseAlarm) + sum(ObjectRejection) + 1));     
                                                
    dp_CombinedStim(isinf(dp_CombinedStim)) = NaN;
    
    ylim = [min([min(dp_CombinedStim,[],'omitnan'),0])-0.1 max([max(dp_CombinedStim,[],'omitnan'),0])+0.1];
                                                    
    clabelstring    = {'Objects'};
    textStartY = max(max(dp_CombinedStim,[],'omitnan')) - 0.15*max(max(dp_CombinedStim,[],'omitnan'));
    textStepY  = max(max(dp_CombinedStim,[],'omitnan'))/20;
    
    ActualdprimeCombinedStimfig = figure('Position',figsize);
    b = bar(dp_CombinedStim);
    b(1).FaceColor = [0.3010 0.7450 0.9330]; % blue
    text(2.0,textStartY,'criterion:','FontSize',18,'fontweight','bold')
    for i = 1:numel(c_CombinedStim)
        text(2.0,textStartY-i*textStepY,clabelstring(i),'FontSize',18)
        if c_CombinedStim(i)<0
    	    text(3.965,textStartY-i*textStepY,num2str(c_CombinedStim(i)),'FontSize',18)
        else
            text(4,textStartY-i*textStepY,num2str(c_CombinedStim(i)),'FontSize',18)
        end
    end
    set(gca,'XTick',1,'XTickLabel',{'Objects'})
    ylabel('dprime (sure only)')
    set(gca,'FontName','Arial','FontSize',24,'LineWidth',2,'Box','off')
    axis([0.5 4.5 ylim])
    title([subjID, ' ', RetrievalDelay ' Actual Discrimination Index - Combined Stim/NoStim (out of sure only responses)'])
    
    fprintf('Saving d prime with combined stim figure...\n')
    savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', RetrievalDelay, SureResponseString,'_Loglinear', '_ActualSuredprimeCombinedStim.png'));
    
    saveas(ActualdprimeCombinedStimfig, [savefile(1:end-4) '.png'],'png');
    fprintf('Done\n')


    
   

end



end


%% Categorize key press responses in the context of image type and stimulation condition

function output = CheckBAConditions(data,novelty,stim,Yes_No_AND_ConfidenceRating)

if strcmp(data{3},novelty) && data{2} == stim && any(ismember(Yes_No_AND_ConfidenceRating,data{9}))
    output = 1;
else
    output = 0;
end

end

function [TotalNoveltyStim] = GetBA_MaxPossible(data,novelty,stim)

for i = 1:size(data,1)
    NoveltyPossible(i) = contains(data{i,3},novelty);
    
    StimPossible(i)   = isequal(data{i,2},stim);
end

NoveltyStimPossible  = NoveltyPossible.*StimPossible;
TotalNoveltyStim     = sum(NoveltyStimPossible);


end
