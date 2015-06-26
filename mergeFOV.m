function funout = mergeFOV(FOVset,expdate)
% this program merges FOVs for data display
%outfile = 'FOV';

% CTper - catch tubes (channels) per field of view. 28 for all experiments so far
% 
% FOVset - fields of view that make up the experimental dataset. Haven't really used it, 	but good information to keep track of.
% 
% ar - Area of fluorescence measurement (roughly equivalent to area of the cell). 		Variable is empty unless using fluorescence data. Not currently used for 		anything downstream.
% 
% d - distance between old pole and new pole
% 
% dvlc - division location. Index location of cell divisions.
% 
% dvln - length of the cell at division (generally shortest length of a division cycle). 		Can get this by plugging dvlc into d.
% 
% dvtm - time of the cell division. Can get by plugging dvlc into t.
% 
% expdate - experimental date
% 
% fli - Average intensity of fluorescence in area ar. Variable is empty unless using 		fluorescence data. Not currently used for anything downstream.
% 
% fstatus - vector showing final status of all cells. 
% 	Index of vector = CTper*(FOV-1)+channel, 
% 	eg: index 29 is for the cell in FOV2, channel 1.
% 	KEY
% 	-1: lost/ejected
% 	 0: empty/never filled
% 	 1: dies
% 	 2: survives
% 
% ndv - number of divisions observed for each cell (regardless of fate). Index is the same 	as for fstatus.
% 
% outfile - name of the file, not needed, but hard to get rid of as a variable.
% 
% t - time point giving the time since the beginning of the experiment.


d = [];                 % raw distance data [pixels]
dvlc = [];              % indices of division events
dvln = [];              % length of cell after each division [pixel]
dvtm = [];              % time of each division [hours]
ndv = [];               % number of divisions of each cell before end of experiment
t = [];                 % raw time data
fli = [];               % fluorescence channel means
ar = [];                % fluorescence areas (px^2)
statind = [];           % index for status used to pick cells in fstatus
for i = FOVset
    %display(i)
    infile = ['FOV_',num2str(i)];
    a = load(infile);
    
    d = [d,a.d];                % raw distance data [pixels]
    dvlc = [dvlc, a.dvlc];          % indices of division events
    dvln = [dvln,a.dvln];     % length of cell after each division [pixels]
    dvtm = [dvtm,a.dvtm];          % time of each division [hours]
    ndv = [ndv,a.ndv];            % number of divisions of each cell before end of experiment
    t = [t,a.t];                % raw time data
    fli = [fli,a.fli];            %fluorescence channel means
    ar = [ar,a.ar];              % fluorescence areas (px^2)
    %outfile = [outfile,'_',num2str(i)];
    CTper = a.CTper;
    clear a
    statind = [statind, (i-1)*CTper+(1:CTper)];
end

% Import data for fstatus, final status of cells
if exist('final_state.mat')==2
    load('final_state.mat');
    fstatus = reshape(fstatus,1,(size(fstatus,1).*size(fstatus,2)));
else
    N = max(FOVset).*CTper; % number of possible statuses
    fstatus = cell2mat(textscan(fopen('final_state.txt'),'%f',N))';
    fstatus = fstatus(statind); % fstatus of cells in merged FOVs only
end

% Purge unwanted variables
clear('N','statind','i','infile')
fstatus((ndv>=nanmax(ndv(fstatus==1)))&(fstatus==-1))=2; % changes "lost" after last death to "survives"
outfile = ['FYLM_',num2str(expdate)];%,'_',outfile];
save(outfile);
funout = outfile;
end