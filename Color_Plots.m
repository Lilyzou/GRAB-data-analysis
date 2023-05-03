close all;
Session_Struct;
frameRate = 15.44;
secResolution = 1;
timeFrame = frameRate * 0.5;
mkdir('Session-Analysis-base');
mkdir('Session-Analysis-500ms');
mkdir('Session-Analysis-frames');
currFolder = pwd;
saveTo1 = strcat(currFolder,'\Session-Analysis-base\');
saveTo2 = strcat(currFolder,'\Session-Analysis-500ms\');
saveTo3 = strcat(currFolder,'\Session-Analysis-frames\');
figNum = 0;
i = 1;
while isempty(session(i).baseData)
    session(i) = [];
end
data = session.baseData;
dataNames = fieldnames(data);
normalizeName = fieldnames(session);
for k = 1:length(fieldnames(session))
    normType = normalizeName{k};
    
    if k == 1
        xlimits = [0,5.55*frameRate];
        xsend = 0.5*frameRate;
        saveTo = saveTo1;
    elseif k == 2
        data = session.timeNormData;
        xlimits = [-0.5*frameRate,2.05*frameRate];
        xsend = 0;
        saveTo = saveTo2;
    else
        data = session.frameNormData;
        xlimits = [0,5.55*frameRate];
        xsend = 0.5*frameRate;
        saveTo = saveTo3;
    end
    ogXlims = xlimits;
    for i = 10:(length(fieldnames(data))-1)
        clear legendCell;
        toskip = 0;
        if (i == 10) || (i == 11)
            xlab = 'time from first reward lick(s)';
            legTitle = 'First Answer Lick';
        elseif (i == 16)
            xlab = 'time from true first lick(s)';
            legTitle = 'True First Lick';
        elseif ((i > 11) && (i < 17))
            xlab = 'time from first lick(s)';
            legTitle = 'First Lick';
        elseif (i == 17) || (i == 18) || (i == 19) || (i == 20)
            xlab = 'time from pole down(s)';
            legTitle = 'Pole Down';
        else
            xlab = 'time from last lick(s)';
            legTitle = 'Last Lick';
        end
        if (i >= 17 && i <= 22) && (k ~= 2)
            xlimits = [ogXlims(1), ogXlims(2)-frameRate/2];
        else
            xlimits = ogXlims;
        end
        
        if (i == 10) || (i == 12) || (i == 16) || (i == 17) || (i == 21)
                trialtype = 'Hit Trials';
        elseif (i == 11) || (i == 13) || (i == 18) || (i == 22)
                trialtype = 'FA Trials';
        elseif (i == 14) || (i == 19)
                trialtype = 'Miss Trials';
        else
                trialtype = 'CR Trials';
        end
        
        figNum = figNum + 1;
        figure(figNum);
        
        
        
        for n = 1:length(session) 
            data = session(n).(normType);
            trace = data.(dataNames{i}){1,1};
            trace = trace(~isnan(trace));
            if (isempty(trace))
                toskip = toskip + 1;
                continue
            end
            minlength(n) = length(trace);
            mintrace(n) = min(trace);
            maxtrace(n) = max(trace);
            traceMean(n,1:min(minlength)) = trace(:,(1:min(minlength)));
            xs = 0:(frameRate/secResolution):(min(minlength)+xsend);
            ax = gca;
            ax.YAxis.Exponent = -3;
            
            ylabel('\DeltaF/F')
            xlim(xlimits);
            xlabel(xlab);
            set(gca,'xtick',xs,'xticklabel',-1:(1/secResolution):length(xs),'FontSize',20,'fontweight','bold');
            if (k == 2)                
                set(gca,'xtick',xs,'xticklabel',0:(1/secResolution):length(xs),'FontSize',20,'fontweight','bold');
            end
            colors = (length(session)-n)/(length(session)+1);
            
            if (i == 10) || (i == 12) || (i == 16) || (i == 17) || (i == 21)
                colorGrad = [colors, colors,1];
            elseif (i == 11) || (i == 13) || (i == 18) || (i == 22)
                colorGrad = [colors, 1, colors];
            elseif (i == 14) || (i == 19)
                colorGrad = [colors, colors, colors];
            else
                colorGrad = [1, colors, colors];
            end
            if k == 2
                xvals = [1:length(trace)];
                xvals = xvals - frameRate/2;
            else
                xvals = [1:length(trace)];
            end
            plot(xvals, trace, 'Color', colorGrad,'linewidth',3);
            hold on;
            legendCell{n-toskip} = comboSave(5*(n-1)+1:4+5*(n-1));
        end
        axis manual;
        if (k == 2)
            stimuli = 0 * ones(1,2);
        else   
            stimuli = frameRate * ones(1,2);
        end
        plot(stimuli,[3*min(mintrace),3*max(maxtrace)],'-.k','linewidth',3);
        hold on;
        if k == 2
            xvals = [1:length(nanmean(traceMean))];
            xvals = xvals-frameRate/2;
        else
            xvals = [1:length(nanmean(traceMean))];
        end
        plot(xvals,nanmean(traceMean),'-k','linewidth',3);
        legendCell{length(session)+1-toskip} = legTitle;
        legendCell{length(session)+2-toskip} = 'Average';
        legend(legendCell,'FontSize',8);
        set(gcf,'Position',[100,100,1200,900]);
        saveas(figNum,strcat(saveTo,strcat(trialtype, " ", legTitle)),'png');
    end
    
    for i = 3:2:9
        clear legendCell;
        toskip = 0;
        xlab = 'time from pole onset(s)';
        legTitle = 'Pole Onset';
        
        if i == 3
            trialtype = 'Hit Trials';
        elseif i == 5
            trialtype = 'Miss Trials';
        elseif i == 7
            trialtype = 'CR Trials';
        else
            trialtype = 'FA Trials';
        end
        
        figNum = figNum + 1;
        figure(figNum);
        for n = 1:length(session)
            
            data = session(n).(normType);
            trace = data.(dataNames{i}){1,1};
            trace = nanmean(trace);
            trace = trace(~isnan(trace));
            
            onset = data.(dataNames{23}){1,1};
            
            if (isempty(trace))
                toskip = toskip + 1;
                continue
            end
            minlength(n) = length(trace);
            mintrace(n) = min(trace);
            maxtrace(n) = max(trace);
            traceMean(n,1:min(minlength)) = trace(:,(1:min(minlength)));
            xs = (onset-frameRate):frameRate/secResolution:min(minlength)+xsend;
            set(gca,'xtick',xs,'xticklabel',-1:(1/secResolution):length(xs),'FontSize',20,'fontweight','bold');
            if (k == 2)                 
                xs = (onset-frameRate):frameRate/secResolution:min(minlength)+xsend;
                set(gca,'xtick',xs,'xticklabel',0:(1/secResolution):length(xs),'FontSize',20,'fontweight','bold');
            end
            xlim(xlimits);
            ax = gca;
            ax.YAxis.Exponent = -3;
            xlabel(xlab)
            ylabel('\DeltaF/F')
            colors = (length(session)-n)/(length(session)+1);
            if i == 3
                colorGrad = [colors, colors,1];
            elseif i == 5
                colorGrad = [colors, colors, colors];
            elseif i == 7
                colorGrad = [1, colors, colors];
            else
                colorGrad = [colors, 1, colors];
            end
            if k == 2
                xvals = [1:length(trace)];
                xvals = xvals - frameRate;
            else
                xvals = [1:length(trace)];
            end
            plot(xvals, trace, 'Color', colorGrad,'linewidth',3);
            hold on;
            legendCell{n-toskip} = comboSave(5*(n-1)+1:4+5*(n-1));
        end
        axis manual;
        onset = data.(dataNames{23}){1,1};
        stimuli = onset * ones(1,2);
        if (k == 2)
            stimuli = (onset-frameRate) * ones(1,2); %not frameRate/2 because onset was not aligned to stimuli like the others
        end

        plot(stimuli,[3*min(mintrace),3*max(maxtrace)],'-.k','linewidth',3);
        hold on;
        if k == 2
            xvals = [1:length(nanmean(traceMean))];
            xvals = xvals-frameRate;
        else
            xvals = [1:length(nanmean(traceMean))];
        end
        plot(xvals,nanmean(traceMean),'-k','linewidth',3);
        legendCell{length(session)+1-toskip} = legTitle;
        legendCell{length(session)+2-toskip} = 'Average';
        legend(legendCell,'FontSize',8);
        set(gcf,'Position',[100,100,1200,900]);
        saveas(figNum,strcat(saveTo,strcat(trialtype, " ", legTitle)),'png');

        
    end
end
close all