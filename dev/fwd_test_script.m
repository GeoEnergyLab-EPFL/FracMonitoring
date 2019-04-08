% forward test script 

%% cleanup first, set global parameters
close all
clearvars
home

% data storage location
datastor = 'gel-nas1'; % 'local' if dataset copied to local drive, 'gel-nas1', or 'enacdrives'

%% choose dataset and load acquisition times
switch datastor
    case 'gel-nas1'
        % data path on gel-nas1
        datapath = pathbyarchitecture('gel-nas1');
    case 'enacdrives'
        datapath = pathbyarchitecture('enac1files');
    case 'local'
        [~, username] = system('whoami');
        datapath = ['/home/' username(1:end-1) '/data/'];
end

% 2019 acquisitions
datayear = 19;
% % injection on slate
datamonth = 03;
dataday = 01;
starttime = '121330';
% test on gabbro
%datamonth = 03;
%dataday = 22;
%starttime = '145034';

% data folder name from experiment date
datafold = [num2str(datayear,'%02d') '-' num2str(datamonth,'%02d') '-' ...
    num2str(dataday,'%02d') '/'];

% extract header info from JSON file
fjson = [datapath datafold num2str(starttime) '.json'];
[jsonhdr,myTransducers,myPlattens,myBlock] = load_header(fjson);

%%

% Choose some rays

% brute force - all rays  
n_c=32;
ch=linspace(1,n_c,n_c)';

myMap=zeros((n_c*n_c),2);
for i=1:n_c
    myMap(1+(i-1)*n_c:i*n_c,1)=i;
    myMap(1+(i-1)*n_c:i*n_c,2)=ch;
end

%
fig_b=plotblockwithplattens(myBlock,myPlattens)

mySR=SourceReceiverPairs(myTransducers,myPlattens,myMap(1:32*32,:));
fig_handle=plotdirectrays(mySR,fig_b)
%%

% take all source -receivers pair except the pairs lying in the same
% platten 
mySR=AllexceptintraPairs(myTransducers,myPlattens);

%%
fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(mySR,fig_b)


%%

% take all sources from platten 'B' and all receivers from platten 'T'
mySR=TwoPlattenPairs(myTransducers,myPlattens,'B','T') ; 
 


%%
% take all S-R for opposite plattens
mySR=AllOppositePairs(myTransducers,myPlattens);



%%  Prepare for diffraction forward

mySR1=TwoPlattenPairs(myTransducers,myPlattens,'B','W') ; 
mySR2=TwoPlattenPairs(myTransducers,myPlattens,'B','E') ; 
mySR3=TwoPlattenPairs(myTransducers,myPlattens,'B','N') ; 
mySR4=TwoPlattenPairs(myTransducers,myPlattens,'B','S') ; 

the_comb_map=CombineMaps(mySR1,mySR2);

the_comb_map2=CombineMaps(mySR3,mySR4);
my_map=union(the_comb_map,the_comb_map2,'rows');

SRdiff=SourceReceiverPairs(myTransducers,myPlattens,my_map);
fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(SRdiff,fig_b);

%----- 
marble = IsotropicSolid(2700,60.*1e9,0.25); % vel in m/s
%%%
ell = Ellipse(90,40,[250./2,250./2,250./2],0*pi/180,45*pi/180,0.*pi/180);
fig_handle=plotEllipse(ell,fig_handle);

%% ----
rock=marble;
x_source=SRdiff.XS_XR(:,1:3);

x_rec=SRdiff.XS_XR(:,4:6);

En=ell.Normal;
Elli1=PointsOnEllipse(ell,720);

tic
results=zeros(SRdiff.n_pairs,2);
for p=1:SRdiff.n_pairs
    
% line containing both source and receiver
Lv = x_source(p,:)-x_rec(p,:);

% intersection of line and ellipse planeEn
t = dot(En,[ell.xc  ell.yc  ell.zc]-x_source(p,:))/dot(En,Lv);
Intersec = x_source(p,:)+t*Lv;

dElli = dists(Intersec',Elli1');

Inear = find(dElli<=min(dElli)*1.35);

% forward propagation: find propagation times for these candidates
% compute the distances between source, diffractor and receiver
dSD = dists(x_source(p,:)',Elli1( Inear,:)');
dRD = dists(x_rec(p,:)',Elli1(Inear,:)');
% compute source-diffractor, and diffractor-receiver direction vectors
directSD = directions(x_source(p,:),Elli1(Inear,:));
%directRD = directions(x_rec(p,:),Elli1(Inear,:));

 dt=(dSD+dRD)/rock.Vp;

[TT, Emintt] = min(dt);
% plot(dt);
% hold on;
results(p,:) = [TT Inear(Emintt)];

end
toc

%%
tic
[results2]=diffractionForward(marble,SRdiff,ell);
toc
plot(results(:,1),results2(:,1),'.')
%% plot things
fig_b=plotblockwithplattens(myBlock,myPlattens)

%
fig_handle=plotdirectrays(SRdiff,fig_b);
fig_handle=plotEllipse(ell,fig_handle);

hold on
plot3(results2(:,2),results2(:,3),results2(:,4),'.g','MarkerSize',30);

%%
% tic 
% dsd=pdist2(x_source,Elli1);
% ddr=pdist2(x_rec,Elli1);
% dt=(dsd+ddr)/rock.Vp;
% [res, ks]=min(dt,[],2);
% toc;

