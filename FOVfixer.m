function funout = FOVfixer(outfold,FOV,varargin)
% Allows user to modify the automatically-selected cell division points

loadpath = [pwd,'/'];

if nargin >2
    tpmx = varargin{1};
    tpidx = nan(1,tpmx); % blank time period index
    tptm = tpidx; % blank time period times
    for j = 1:tpmx
        tfilename = ['tp',num2str(j),'-fov',num2str(FOV-1), '.txt'];
        %tstmp = cell2mat(textscan(fopen(filename),'%f%f'));
        tstmp = load([loadpath,'timestamp/',tfilename]);
        tpidx(j) = tstmp(end,1);
        tptm(j) = tstmp(end,2);
    end
end

fate = load([loadpath,'summary/final_state.txt']);
celfat = fate(FOV,:);
inoutfile = [loadpath,outfold,'/','FOV_',num2str(FOV)];
load(inoutfile);
dvrw = size(dvtm,1);%round(max(max(t))./rfin); % rows of matrix for recording divisions
i = 1;
while ((i>0)&&(i<29))
%plot(t(:,i),d(:,i),'ok',dvtm(:,i),dvln(:,i),'r.',t(:,i),fli(:,2.*(i-1)+1).*500,'y',t(:,i),fli(:,2.*(i)).*500,'m')
%plot(t(:,i),d(:,i),'ok',dvtm(:,i),dvln(:,i),'r.')
plot(t(:,i),d(:,i),'ok')
hold on
plot(dvtm(:,i),dvln(:,i),'r.','MarkerSize',14)
hold off
if nargin >2
    hold on
    stem(tptm./3600,50.*ones(tpmx,1),'g')
    hold off
end
boxx = xlim;
boxx = boxx(2);
boxy = ylim;
boxy = boxy(2);
boxw = boxx./10;
boxh = boxy./10;
% Next or Previous
rectangle('Position',[boxx-boxw,boxy-boxh,boxw,boxh])
text(boxx-boxw/1.5,boxy-boxh/2,'Next')
rectangle('Position',[boxx-boxw-boxw,boxy-boxh,boxw,boxh])
text(boxx-boxw-boxw/1.5,boxy-boxh/2,'Prev')
text(boxx-(4*boxw),boxy-boxh/2,['Channel ',num2str(i)])
% Add or Delete
rectangle('Position',[boxx-boxw,boxy-boxh-boxh,boxw,boxh])
text(boxx-boxw/1.5,boxy-boxh-boxh/2,'Add')
rectangle('Position',[boxx-boxw-boxw,boxy-boxh-boxh,boxw,boxh])
text(boxx-boxw-boxw/1.5,boxy-boxh-boxh/2,'Del')
rectangle('Position',[boxx-boxw,boxy-boxh-boxh-boxh,boxw,boxh])
text(boxx-boxw/1.1,boxy-boxh-boxh-boxh/2,['Fate: ',num2str(celfat(i))])
%text(boxx-(4*boxw),boxy-boxh-boxh/2,['Fate: ',num2str(celfat(i))])

[PNx,PNy] = ginput(1);
if ((PNx>(boxx-boxw))&&(PNy>(boxy-boxh))) % Advances to the next cell
    i = i+1;
elseif ((PNx>(boxx-boxw))&&(PNy<(boxy-boxh))&&(PNy>(boxy-boxh-boxh))) %ADD
    dloc1 = dvlc(:,i);
    text(boxx-boxw-boxw-boxw-boxw,boxy-boxh-boxh-boxh,'click to add, "Enter" to continue')
    apt = ginput; %points to be added
    for rr = 1:size(apt,1) % This loop adds the points clicked
       adptt = (abs(t(:,i)-apt(rr,1))); % finds difference between point to be added and all time points
       adptd = (abs(d(:,i)-apt(rr,2))); % finds difference between point to be added and all length points
       adpt = adptt+adptd;
       adti = find(adpt == min(adpt)); % finds smallest difference
       dloc1 = [dloc1; adti];
    end
    dloc1 = sort(dloc1); % sorts into order
    dloc1 = dloc1(dloc1>0); % gets rid of NaN
    ndv(i) = length(dloc1); %%% updated number of divisions
    dvtm(:,i) = nan(dvrw,1); % blank current channel division times
    dvln(:,i) = nan(dvrw,1); % blank current channel lengths
    dvlc(:,i) = nan(dvrw,1); % blank current channel division indices
    dvtm(1:ndv(i),i) = t(dloc1,i); % updated division times
    dvln(1:ndv(i),i) = d(dloc1,i); % updated division lengths
    dvlc(1:ndv(i),i) = dloc1; % updated division indices
%    beep
elseif ((PNx<(boxx-boxw))&&(PNy<(boxy-boxh))&&(PNx>(boxx-boxw-boxw))) %REMOVE
    dloc1 = dvlc(:,i);
    text(boxx-boxw-boxw-boxw-boxw,boxy-boxh-boxh-boxh,'click to remove, "Enter" to continue')
    dpt = ginput; %points to be removed
    for rr = 1:size(dpt,1) % This loop removes the points clicked
       dlpt = (abs(dvtm(:,i)-dpt(rr,1))); % finds difference between point to be removed and all points
       dlti = find(dlpt == min(dlpt)); % finds smallest difference
       dloc1(dlti) = NaN; % replaces value with NaN
    end
    dloc1 = dloc1(dloc1>0); % removes NaN
    ndv(i) = length(dloc1); %%% updated number of divisions
    dvtm(:,i) = nan(dvrw,1); % blank current channel division times
    dvln(:,i) = nan(dvrw,1); % blank current channel lengths
    dvlc(:,i) = nan(dvrw,1); % blank current channel division indices
    dvtm(1:ndv(i),i) = t(dloc1,i); % updated division times
    dvln(1:ndv(i),i) = d(dloc1,i); % updated division lengths
    dvlc(1:ndv(i),i) = dloc1; % updated division indices
elseif ((PNx>(boxx-boxw))&&(PNy<(boxy-boxh-boxh))&&(PNy>(boxy-boxh-boxh-boxh))) %FATE
    if celfat(i)==2
        celfat(i) = -1;
    else
        celfat(i) = celfat(i)+1;
    end
else
    i = i-1;
end % if loop

save(inoutfile);
funout = inoutfile;%load(outfile);

end % while loop

% save(outfile);
% funout = outfile;%load(outfile);

end % function

