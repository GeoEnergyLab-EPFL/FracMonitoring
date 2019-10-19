function [energy, SRmap]  = energy_evolution(dataseq1,endnoise,myTransducers,myPlattens,varargin)
% Dong Liu -- 17/09/2019
% function creating energy changes compared with the reference signal for
% all transmisson transducer pairs on the top and bottom platten or for the
% source or receiver and its corresponding receivers or sources on the
% opposite platten
% optional argument 1: number of the transducer, corresponding to the SRmap
% optional argument 2: type of the transducer, when having this one plot
% automatically the one transducer S or R-multi corresponding transducers R
% or S
% when option1 and option 2 are both[], one deals automatically with the
% T-B all facing pairs
% optional argument 3,4,5,6: calculate partially the wave energy(norm)
% reference sequence, [Vp Vs], Fs, actual half window size
nargin = length(varargin);

if ~isempty(varargin) && ~isempty(varargin{1}) && ~isempty(varargin{2})
    [energy, SRmap] = evolution_single(dataseq1,endnoise,myTransducers,myPlattens,varargin{1},varargin{2},[],[],[],[]);
    if nargin>=6 && ~isempty(varargin{4}) && ~isempty(varargin{5}) && ~isempty(varargin{6})
        [energy, SRmap] = evolution_single(dataseq1,endnoise,myTransducers,myPlattens,varargin{1},varargin{2},...
            varargin{3},varargin{4},varargin{5},varargin{6});
    end
else
    [energy, SRmap] = evolution_transmit(dataseq1,endnoise,myTransducers,myPlattens,[],[],[],[]);
    if nargin>=6 && ~isempty(varargin{4}) && ~isempty(varargin{5}) && ~isempty(varargin{6})
        [energy, SRmap] = evolution_transmit(dataseq1,endnoise,myTransducers,myPlattens,...
            varargin{3},varargin{4},varargin{5},varargin{6});
    end
end

end

function [energy_transmit, SRmap] = evolution_transmit(dataseq1,endnoise,myTransducers,myPlattens,refseqnb,VArray,Fs,windowsize)

% detect all the transducers on the top and the bottom, regardless of the
% srouce or receiver
sidemarker='T'; % we plot the top platten by default 
for p=1:length(myPlattens)
    if myPlattens(p).face==sidemarker
        id_p=p;
    end
end
id_transducer=find(myTransducers.platten == myPlattens(id_p).id);
transmitnumber=myTransducers.channel(id_transducer)+1;
% build the transducer Pairs and calculate the energy change
transmitPairs=SourceReceiverPairs(myTransducers,myPlattens,[transmitnumber transmitnumber]);
energy_transmit  = energy_calculation(dataseq1,endnoise,transmitPairs,refseqnb,VArray,Fs,windowsize);
SRmap = transmitPairs.SRmap;
end


function [energy_all, SRmap] = evolution_single(dataseq1, endnoise, myTransducers,myPlattens,number,type,refseqnb,VArray,Fs,windowsize)
% now only working for the top and bottom platten
% we locate the transducer position with the number and the type
id_transducer = find((myTransducers.channel+1==number).*(myTransducers.type==type));
for p=1:length(myPlattens)
    if myTransducers.platten(id_transducer)==myPlattens(p).id
        oposide=myPlattens(p).face;
    end
end
switch oposide
    case 'N'
        sidemarker='S';
    case 'S'
        sidemarker='N';
    case 'E'
        sidemarker='W';
    case 'W'
        sidemarker='E';
    case 'T'
        sidemarker='B';
    case 'B'
        sidemarker='T';
end 

if type=='S'
    side_source=oposide;
    side_receiver=sidemarker;
else
    side_source=sidemarker;
    side_receiver=oposide;
end
oppositeSRs = SinglePlattenPairs(myTransducers,myPlattens,side_source,side_receiver,number,type);
energy_all  = energy_calculation(dataseq1,endnoise,oppositeSRs,refseqnb,VArray,Fs,windowsize);
SRmap = oppositeSRs.SRmap;

end