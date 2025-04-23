close all
clear all
SonjasDaten=0 % put 1 if data from Sonja is used
% Do the following in Igor
% Go to root:saved_traces:FRET:offset:FRET:corrected
% then put in command line:
% string list
% list=WaveList("*_x*_y*", ";", "" )
% save/B/J/W list as "OutputFile.txt"
%
% This codes does the follwing
% 1. import txt with all traces (they are pairs of 3: r_r, r_g, g_g)
%     r_r = acceptor fluorescence after acceptor excitation
%     r_g = acceptor fluorescence after donor excitation
%     g_g = donor fluorescence after donor excitation
% 2. calculate FRET efficiency from g_g and r_g signal 
% 3. individual FRET traces are saved to a folder
%
% What you have to do: 
% 1. input the folder (uncomment line 24) and filename (line 26) where your data file
% from Igor is located in the following lines
%% import data
if SonjasDaten==0
    %cd Z:\path\to\folder
    % % go to folder with your big txt-file and put its name here:
    filename ='example_data.txt';
    % Auf Reihenfolge der Farben achten!!!
    delimiterIn = '\t';
    headerlinesIn = 1;
    A = importdata(filename,delimiterIn,headerlinesIn);

    % select values
    Data=A.data(:,1:end); % for all traces, meaning r_r r_g g_g
    Data(:,1:3:end)=[]; % for g_g and r_g run line above and this line
    colorsselected=2;
    tracenames=A.textdata(1,2:3:end);%
    T=[0:0.5:200.5]';% time in seconds
    T=T(1:size(Data,1));% to make time string the same length as the longest trace
end
%% import Sonjas data 
if SonjasDaten==1
    folder='Z:\_personalDATA\JS+LV_4F-TIRF\Dynamics\UdoSeifert\FRETtraces_forUdoSeifert\elife-57180-fig2-data1-v1\smFRETdatasets\crowding\Ficoll20';
    filename='Sonja_crowding_20percentFicoll____';%is not the true filename, but the savename
    files= dir(fullfile(folder,'*.dat'));
    mydir  = folder;
    idcs   = strfind(mydir,'\');
    % Auf Reihenfolge der Farben achten!!!
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 3);
    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = "\t";
    % Specify column names and types
    opts.VariableNames = ["d_d", "a_d", "a_a"];
    opts.VariableTypes = ["double", "double", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    Data= nan(750,length(files)*3);
    % Import the data
    for i = 1:length(files)
        j=i*3;
        % Import the data
        trace_tab = readtable(fullfile(files(i).folder,files(i).name),opts);
        Data(1:height(trace_tab),j-2:j) = table2array(trace_tab); % Convert to matrix
    end
    % select values sonjas Reihenfolge <d_d>	<a_d>	<a_a>
    Data(:,3:3:end)=[]; % for g_g and r_g run  this line
    colorsselected=2;
    tracenames={files.name};%
    T=[0:0.2:150-0.2]';% time in seconds
end
%% calculate FRET E
calcFRET=1;
if calcFRET && SonjasDaten==0
    FRETE=[];
    for i=1:2:(size(Data,2)-1) % calculate FRET E
        FRETE(:,i)=Data(:,i)./(Data(:,i)+Data(:,i+1));
    end
    disp('blue')
elseif calcFRET && SonjasDaten==1
    FRETE=[];
    for i=1:2:(size(Data,2)-1) % calculate FRET E
        FRETE(:,i)=Data(:,i+1)./(Data(:,i)+Data(:,i+1));
    end
    disp('green')
else
        FRETE=[];
    for i=1:2:(size(Data,2)-1) % calculate FRET E
        FRETE(:,i)=Data(:,i+1);
    end
    disp('red')
end
%delete Zeros
FRETE=FRETE(:,1:2:end);
%% saving traces
saveto= strcat(filename(1:end-4),'_individual_FRETtraces');
mkdir(saveto)

% run a loop to save each column into text file
for i = 1:size(Data,2)/colorsselected; 
    
    %     k=reshape( repmat( [1:size(Data,2)/colorsselected], colorsselected,1 ), 1, [] ); % to get 1 1 1 2 2 2 ....
    filenamearray = strcat(tracenames(i),'.txt') ; %strcat(num2str(k(i)),'_',tracenames(i),'.txt') ; % filename of saved trace
    filename=cat(2,filenamearray{:});
    if SonjasDaten
        filename=strrep(filename,'.dat','');
        filename=strrep(filename,'o_g','FRETE');
    end
    filename=strrep(filename,'r_g','FRETE');
    data_col = FRETE(:,i) ; %select data of trace
   
    data_col=rmmissing(data_col); 
    save(strcat(saveto,'/',filename),'data_col','-ascii') ;
end
%% plot exemplary traces

traceNumber=13;%uneven number!
num_no_nans = sum(  ~isnan(Data(:,traceNumber)) );
figure
%size
x0=10;
y0=10;
width=24.6;
height=8;
set(gcf,'units','centimeters','position',[x0,y0,width,height])
hold on
% LEFT
yyaxis left
plot(T,Data(:,traceNumber),'r')%plot(Data(:,traceNumber),'r')
plot(T,Data(:,traceNumber+1),'color','green','LineStyle','-','LineWidth',1)%plot(Data(:,traceNumber+1),'color','green','LineStyle','-','LineWidth',1)
ylabel('fl. intensity / a.u.')
ylim([-9500 25000])
% RIGHT
yyaxis right
plot(T,FRETE(:,(traceNumber+1)/2),'k')%plot(FRETE(:,(traceNumber+1)/2),'k')
ylim([-0.2 1.2])
ylabel('FRET E')
ax = gca;
% formatting
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% set(ax.YAxis(2), 'TickValues', [0 [] 1 []])
newYTick=ax.YAxis(1).TickValues/1000;
set(ax.YAxis(1), 'TickLabels', newYTick)
% legend('Aem Dex','Savitzky-Golay smooth with window=5','mid line','ground truth','state allocation','Box','on')
xlabel('seconds')

legend('A_{em}D_{ex}','A_{em}D_{ex}','FRET efficiency','Box','off')
% title(strcat('trace',' ',num2str(traceNumber)))
% fontsize(14,'points')
set(gca, 'XTickLabels', ax.XTick/2)
% xlim([0 82*2])
%%
figure
histogram(FRETE,[-5:0.04:5])
xlim([-0.5,1.5])
xlim([-1,1.5])