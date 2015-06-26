%%% Function for crunching one FOV of raw data from FYLM_critic 2.0
%%%
%%% assumes current folder is folder with input text files. Text files
%%% should have two columns for no fluorescence and an additional five
%%% columns per fluorescence channel
%%%
%%% FYLM Critic output data:
%%% separate files for each catch channel
%%% no separate files for Time Periods
%%% times are all relative to beginning of experiment
%%% fluorescence data is alphabetical by name in Elements

function outstruct = FYLMcrunchStart(FOV,fluorescence)
% if ~isstr(expdate)
%     error('expdate must be a text string')
% end
% if ~isstr(path)
%     error('path must be a text string')
% end



CTper = 28; % Total number of catch channels per FOV

fmt = '%f%f';
if fluorescence>0
    fct = 1:fluorescence:(CTper.*fluorescence);
    for fi = 1:fluorescence
        fmt = [fmt,'%f%f%f%f%f'];
    end
else
    fli = [];
    ar = [];
end

for ct = 1:CTper
    
    filename = [num2str(FOV), '_', num2str(ct), '.txt'];
    fylm = cell2mat(textscan(fopen(filename),fmt));

    t(:,ct) = fylm(:,1);    % time of each frame
    d(:,ct) = fylm(:,2);    % pixel length of cell in each frame
    if fluorescence>0       
        fli(:,fct(ct):(fct(ct)-1+fluorescence)) = fylm(:,3:5:end); % normalized fluorescence intensity in each frame, by fluorophore
        ar(:,ct) = fylm(:,8);
    end
end

t = t/3600; % converts to hours;
mpkd = 40; % minimum peak distance index separation (~80 minutes)
rfin = (mpkd.*2)./60; % minimum division interval time in hours (80 minutes)
dvrw = round(max(max(t))./rfin); % rows of matrix for recording divisions
dvlc = nan(dvrw,CTper); %preallocate matrix for division locations
dvtm = nan(dvrw,CTper); %preallocate matrix for division times
dvln = nan(dvrw,CTper); %preallocate matrix for division lengths
ndv = nan(1,CTper); % preallocated matrix for number of divisions
cterr = 0; % count of quantification errors resulting in unrealistically small lengths
for di = 1:CTper % Cleanup of raw data
    for dj = 2:length(d)
        if (d(dj,di).*5)<(d(dj-1,di))
            cterr = cterr+1;
            d(dj,di) = NaN; % Gets rid of unrealistically small datapoints caused by quantification errors
            if fluorescence>0
                fli(dj,fct(di):(fct(di)-1+fluorescence)) = NaN;
            end

        end
        if isnan(d(dj,di))
            d(dj,di)= d(dj-1,di); % fills in NaN values in distance with the last known distance for all catch channel
            if fluorescence>0
                fli(dj,di)= fli(dj-1,di);
            end

        end
    end
    dsmooth = medfilt1(d(:,di),9);
    [~,dlocs]= findpeaks(diff(dsmooth.*-1),'minpeakdistance',mpkd,'threshold',5); % find division times
    dlocs = dlocs+1;
    dtimes = t(dlocs,di);
    ndv(di) = length(dtimes); %%% This is the number of divisions per cell (RLS for dead cells)
    dvlc(1:ndv(di),di) = dlocs;
    dvtm(1:ndv(di),di) = dtimes;
    dvln(1:ndv(di),di) = d(dlocs,di);
end

outfile = ['FOV_',num2str(FOV)];
save(outfile);
outstruct = outfile;%load(outfile);

%outstruct = {CTper; filename; t; d; fli; cterr};
end

% t = cell2mat(outstruct(3));