% First run the code "calcFRETfromIgor.m" to generate individual FRET
% traces
% This code does the following:
% 1. Import individual FRET traces
% 2. Plot 3D histogram over position in trace and FRET E
% 3. Plot time course
%
% What you have to do:
% 1. out in the file path to the folder where your individual FRET traces
% are stored
% 2. set time axis including time break if needed

clear all
close all
SonjasDaten=0 % put 1 if you used data from Sonja
saveON=1; % put 1 if you want to save your plot
%% set time axis PUT TIME BREAK here
if SonjasDaten
    time=[0:0.2:150-0.2];
else
time=[0:0.5:100-0.5];
end
% set time break as NaN
time=[0:0.5:10,NaN,31:0.5:120-0.5]; % for TIME Break
time=[0:0.5:10,NaN,15:0.5:103.5]; % for 5s TIME Break after 21 loops
%% Import data FRET traces
folder ='example_data_individual_FRETtraces';
%folder='Z:\path\to\datafolder';

files= dir(fullfile(folder,'*.txt'));
 mydir  = folder;
 idcs   = strfind(mydir,'\');
%  folderup = extractBefore(mydir,idcs(end));
%  saveto =fullfile(folderup,'traces4Fabianideal'); %name of folder where traces should go

%  if ~exist(saveto, 'dir')
%    mkdir(saveto)
% end
stringparts=split(folder,"\") % helps with naming when saving plot
titlename=stringparts{end};
titlename2=titlename(1:end-22);
titlename=replace(titlename2,'_',' ');
%% prepare readtable
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = "e01";
opts.VariableTypes = "double";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";
% create empty trace array
if SonjasDaten
traces= nan(750,length(files));
else
    traces= nan(200,length(files));
end
% actual trace reading and saving loop
for i = 1:length(files)
    % Import the data
    trace_tab = readtable(fullfile(files(i).folder,files(i).name),opts);
    trace_array=table2array(trace_tab);
    traces(1:height(trace_tab),i) = trace_array; % Convert to matrix
    % first row of traceinfo stores length of trace
    traceinfo(1,i)=length(trace_array); 
    % check if trace bleaches in high or low FRET
    if length(trace_array)>3
        if mean(trace_array(end-2:end))>0.5
            % second row of traceinfo stores from where trace bleaces
            traceinfo(2,i)=1;
        end
    else
        traceinfo(2,i)=nan;
    end
end
%% histograms per time point
edgesinput=[-2:0.04:2]; % histogram bin edges
edgesbarchart=[(-2+0.04/2):0.04:(2-0.04/2)]; % histogram bin edges when plottet as bar chart
for i=1:size(traces,1)
   [N(i,:),edges] = histcounts(traces(i,:),edgesinput) ;
end
Nnorm=N./sum(N,2); % normalize
Nnormsum=sum(Nnorm,2);%to test if always 1
Nsum=sum(N,2);% to test how many traces at which point in time

% Histogram of all values, including outliers, because when calculating FRET E g_g-Intensity is sometimes divided by 0
edgesinputgross=[-500:0.04:500];
for i=1:size(traces,1)
   [Ngross(i,:),edgesgross] = histcounts(traces(i,:),edgesinputgross) ;
end
Ngrossnorm=Ngross./sum(Ngross,2);
Ngrosssum=sum(Ngross,2);
%done separately for traces which end in low/high FRET
highFRETindex=find(traceinfo(2,:)); % find nonzero values
lowFRETindex=find(traceinfo(2,:)==0); % find zero values
tracesHighFRET=traces(:,highFRETindex);
tracesLowFRET=traces(:,lowFRETindex);
for i=1:size(tracesHighFRET,1)
   [NtracesHighFRET(i,:),edgesTrash] = histcounts(tracesHighFRET(i,:),edgesinputgross) ;
end
NtracesHighFRETsum=sum(NtracesHighFRET,2);
for i=1:size(tracesLowFRET,1)
   [NtracesLowFRET(i,:),edgesTrash] = histcounts(tracesLowFRET(i,:),edgesinputgross) ;
end
NtracesLowFRETsum=sum(NtracesLowFRET,2);

% calculate ratio of low FRET
lowFRETRatio=NtracesLowFRETsum./Ngrosssum;

if SonjasDaten
    Npart=Nnorm(1:5:200,:);
else
Npart=Nnorm(1:2:80,:); % only take first 80 frames
end
% plot 3D histogram
figure
bar3(edgesbarchart,Npart',1);
ylim([-1,2])
xlim([0.5,40.5])
zlim([0 0.14])
xlabel('position in trace [s]')
ylabel('FRET efficiency')
zlabel('normalised abundance')
title(titlename)

% set(gca,'XTick',0:10:40)
% set(gca,'XTickLabel',0:5:12)

% "normal edges": Cut of negative and positive outliers which come from
% division by 0 
% "nomal edges" = -1 till 2 of FRET E
low=Nnorm(:,25:62); % -1 to 0.44 positions correspond to positions in edgesinput
high=Nnorm(:,63:100); % 0.52 to 2
lowSum=sum(low,2);
highSum=sum(high,2);
%% plot time course of FRET efficiency
figure % Histogramm over "normal edges"
hold on
% LEFT
yyaxis left
plot(time,lowSum,'color','green')
plot(time,highSum,'color','red','Linestyle','-')
xlabel('position in trace [s]')
ylabel('percentage in respective FRET E bin')
xlim([0 60])
% RIGHT
yyaxis right
% plot(time,Nsum,'b')% damit wirds noisy durch die edges die zeug abschneiden
plot(time,Ngrosssum,'b')
plot(time,NtracesLowFRETsum,'Color','#2d712c','LineStyle','--')
plot(time,NtracesHighFRETsum,'Color','#800000','LineStyle','-.')

ylabel('amount of traces')
% formatting
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
legend('low FRET bin (-1 -> 0.5)','high FRET bin (0.5 -> 2)',...
    'amount of traces',...
    'which bleach from low FRET',...
    'which bleach from high FRET','box','off');
title(titlename)
% save plot as png
if saveON
saveas(gcf,strcat('FREThistosovertime\',titlename2,'_TimeCourse_withBleachingNumbers.png'));
% saveas(gcf,strcat('FREThistosovertime\Aha1data\',titlename2,'_TimeCourse_withBleachingNumbers.png'));
end
%% histogramm over everything, including outliers
lowgross=Ngrossnorm(:,1:12513);
highgross=Ngrossnorm(:,12514:end);
lowgrossSum=sum(lowgross,2);
highgrossSum=sum(highgross,2);
figure %histogramm über alles, auch ausreißer
hold on
plot(time,lowgrossSum,'color','green')
plot(time,highgrossSum,'color','red')
xlabel('position in trace [s]')
ylabel('percentage in respective FRET E bin')
xlim([0 40])
legend('low FRET bin (-{\infty} -> 0.5)','high FRET bin (0.5 -> {\infty} )','box','off')
title(titlename)
%% bleaching time histograms
figure
histogram(traceinfo(1,:),'BinEdges',[0:1:200])
title('length of traces');

figure
histogram(traceinfo(1,highFRETindex),'BinEdges',[0:1:200])
title('length of traces which bleach from high FRET');

figure
histogram(traceinfo(1,lowFRETindex),'BinEdges',[0:1:200])
title('length of traces which bleach from low FRET');

