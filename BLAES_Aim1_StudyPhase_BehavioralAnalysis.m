%Krista Wahlstrom, 2023
%This is the behavioral analysis script for BLAES Aim 1 Long-Delay for the
%Study Phase (encoding) session.


%% StimulusCodes
% 1 - 2314 Images
% 2501-2504 Instructions screens
% 2601 Sync pulse
% 2701 No Stim
% 2702 Stim
% 5200 ISI (w/fixation cross onscreen)
% 2801 - 5115 Variable timed fixation cross
% 5301 - 6319 Confidence ratings (only applies to the TEST PHASE)


%%

clear;
close all;


addpath(genpath(fullfile(cd,'BCI2000Tools')))

subjID         = 'BJH029';


d = dir(fullfile(cd,'data',subjID,'Study','*.dat'));


%Set this to true if patient is responding as you would expect during the
%presentation of the images. Set this to false, if patient did not respond
%during the presentation of the image (and instead made responses during
%the stimulation period 2701/2702 or the inter-trial ISI period 5200 or 
% during the fixation cross period 2801-5115) 
respondDuringImage = true;

%Create response string for saving separate figure files and .mat files based on which analysis is run 
if respondDuringImage == true
    ResponseString = '_RespondDurImage';
else
    ResponseString = '_RespondAftImage';
end


%Combine all .dat files for the study phase session into one matrix
iter2 = 1;
    for file = 1:size(d,1)
        
        [~, states, parameters] = load_bcidat(fullfile(d(file).folder,d(file).name));
        pause(1);

        SEQ = parameters.Sequence.NumericValue;
        KD = states.KeyDown;
        StimCode = states.StimulusCode;

        SEQ(SEQ < 1) = [];
        SEQ(SEQ > 2314) = [];

        %Replaces the stimcode values that are not equal to a SEQ/image
        %stimulus code, with the previously shown image's stimulus code/SEQ
        %value if respondDuringImage is set to "false" so that keypresses
        %can be gathered from during the intertrial period
        if(~respondDuringImage)
            lastValidImage = -1;
            for i = 1:length(StimCode)
                if(lastValidImage ~= StimCode(i) && ismember(StimCode(i), SEQ))
                    lastValidImage = StimCode(i);
                elseif(lastValidImage ~= -1 && ~ismember(StimCode(i), SEQ))
                    StimCode(i) = lastValidImage;
                end
            end
        end

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

        %Create collectStudyData matrix of filename, stimulus code/image
        %code, key press
        for i = 1:length(SEQ)
            collectStudyData{iter2,1} = parameters.Stimuli.Value{2,SEQ(i)};
            collectStudyData{iter2,2} = SEQ(i);
            
            
            idx = ismember(keyPresses(:,2), SEQ(i));
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
                collectStudyData{iter2,3} = 'Indoor';
            elseif pressForImage == 39
                collectStudyData{iter2,3} = 'Outdoor';
            else
                collectStudyData{iter2,3} = 'Non-Response Key';
            end

            collectStudyData{iter2, 4} = pressForImage;
            iter2 = iter2 + 1;
        end

    end
    
    save(fullfile(cd,'data',subjID,'Study',strcat(subjID,'_Study_Data',ResponseString,'.mat')),'collectStudyData')


%% Behavioral Analysis
%Create a copy of collectStudyData so the original data is preserved. 
Indoor_Outdoor_StudyData = collectStudyData;


%Change the Indoor_Outdoor_StudyData matrix into a table
Indoor_OutdoorTable = cell2table(Indoor_Outdoor_StudyData);


%Add all possible behavioral response options (Indoor, Outdoor, and Non-Response Key
% to column 3 of Indoor_OutdoorTable so that in the case where a patient doesn't use
% one of those responses, they're still added to the Counts table below
% to be analyzed.
Indoor_OutdoorTable(size(Indoor_OutdoorTable,1)+1,3) = {'Indoor'};
Indoor_OutdoorTable(size(Indoor_OutdoorTable,1)+1,3) = {'Outdoor'};
Indoor_OutdoorTable(size(Indoor_OutdoorTable,1)+1,3) = {'Non-Response Key'};



%Get the group counts for indoor/outdoor/non responses
%IncludeEmptyGroups will also display values in the table that are zero
Indoor_Outdoor_Counts = groupcounts(Indoor_OutdoorTable,{'Indoor_Outdoor_StudyData3'}, 'IncludeEmptyGroups', true);


%Subtract 1 from each of the groupcounts in the Indoor_Outdoor_StudyData table that's
%associated with a 'outdoor', 'indoor', and 'non-response' value, because we artificially added these
%responses in line 130 above
for l = 1:size(Indoor_Outdoor_Counts,1)
    if contains(Indoor_Outdoor_Counts{l,1}, 'Indoor')
        Indoor_Outdoor_Counts{l,2} = Indoor_Outdoor_Counts{l,2}-1;
    elseif contains(Indoor_Outdoor_Counts{l,1}, 'Outdoor')
        Indoor_Outdoor_Counts{l,2} = Indoor_Outdoor_Counts{l,2}-1;
    elseif contains(Indoor_Outdoor_Counts{l,1}, 'Non-Response Key')
        Indoor_Outdoor_Counts{l,2} = Indoor_Outdoor_Counts{l,2}-1;
    end
end


%Sort IndoorOutdoorCounts alphabetically
Indoor_Outdoor_Counts = sortrows(Indoor_Outdoor_Counts,1);

%Remove the percent column from IndoorOutdoorCounts because it's not needed for
%analysis
Indoor_Outdoor_Counts.Percent = [];



%Create text for the upper corner of the bar plot that lists each
%condition and the number of responses for each
plottext = {strcat('Indoor: ',num2str(Indoor_Outdoor_Counts{1,2})),strcat('Outdoor: ',num2str(Indoor_Outdoor_Counts{3,2})),...
    strcat('Non-Response Key: ',num2str(Indoor_Outdoor_Counts{2,2}))};




%Plot each behavioral response condition as a percentage of the total
%images in that condition (i.e. percentage of items/scenes/scrammbled)

b = bar([(Indoor_Outdoor_Counts{1,2}/(Indoor_Outdoor_Counts{1,2}+Indoor_Outdoor_Counts{2,2}+Indoor_Outdoor_Counts{3,2}))*100;...
    (Indoor_Outdoor_Counts{3,2}/(Indoor_Outdoor_Counts{3,2}+Indoor_Outdoor_Counts{1,2}+Indoor_Outdoor_Counts{2,2}))*100;...
    (Indoor_Outdoor_Counts{2,2}/(Indoor_Outdoor_Counts{2,2}+Indoor_Outdoor_Counts{1,2}+Indoor_Outdoor_Counts{3,2}))*100]);

xlim([0 5])
ylim([0 105])

b.FaceColor = 'flat';
%Set colors for Indoor, Outdoor, Non-response bars
b.CData(1,:) = [1 0 0]; %red
b.CData(2,:) = [0 0 1]; %blue
b.CData(3,:) = [0 0 0]; %black


%Set the plottext labels to be positioned above bar 3.6 on the xaxis, and the
%yaxis position to be at 70% of the max number of behavioral responses that
%occur
text(3.6,max(Indoor_Outdoor_Counts{:,2})*0.70, plottext)
title([subjID, ' ', 'Study Phase Responses'],'fontweight','bold','fontsize',16)
xticklabels({'Indoor', 'Outdoor','Non-Response'})
xlabel('Behavioral Response','fontweight','bold','fontsize',12)
ylabel('Counts (% of images shown)','fontweight','bold','fontsize',12)

%get current figure "gcf" for saving purposes
f = gcf;

%Save bar graph
savefile = fullfile(cd, 'figures', subjID, strcat(subjID, '_', ResponseString, '_StudyPhase_BehavioralResponses.png'));
saveas(f, savefile);