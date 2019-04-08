
function [ptop, pbottom, pnorth, psouth, peast, pwest, ...
    dcross, dhoriz, dvert1, dvert2, ddiag1, ddiag2,...
    Xsources, Ysources, Zsources, Sources, Xreceivers, Yreceivers, Zreceivers, Receivers,...
    SStmp, RRtmp, SS, RR]...
    = forward_ellipse(debug, Elli1)

% define source and receiver positions
% platen position wrt center of block
ptop = 125;
pbottom = -125;
pnorth = 125;
psouth = -125;
peast = 125;
pwest = -125;

% spacing on platen
dcross = 18;    % top platen cross spacing
dhoriz = 50;    % side platen horiz spacing
dvert1 = 20;    % side platen first vert spacing
dvert2 = 30;    % side platen farther vert spacing
ddiag1 = 34;    % top platen first diag spacing
ddiag2 = 28;    % top platen additional diag spacing

% position of source transducers (at center) %%% WHAT?
Xsources = [zeros(1,8), [4 3 2 1]*dcross, -[1 2 3 4]*dcross,...
    ones(1,6)*peast, -dhoriz , 0, dhoriz, dhoriz, -dhoriz, 0,...
    [ddiag1+ddiag2 ddiag1 -ddiag1 -ddiag1-ddiag2]/sqrt(2)];
Ysources = [[4 3 2 1]*dcross, -[1 2 3 4]*dcross, zeros(1,8),...
     dhoriz , 0, -dhoriz, -dhoriz, 0, dhoriz, ones(1,6)*pnorth,...
     [ddiag1+ddiag2 ddiag1 ddiag1 ddiag1+ddiag2]/sqrt(2)];
Zsources = [ones(1,16)*ptop, ones(1,3)*-dvert1, ones(1,3)*dvert1,...
    ones(1,3)*-(dvert1+dvert2), ones(1,3)*(dvert1+dvert2),...
    ones(1,4)*ptop];
% combine all three coordinates
Sources = [Xsources; Ysources; Zsources];

% position of receiver transducers (at center)
Xreceivers = [Xsources(1:16), ones(1,6)*pwest, Xsources(23:end)];
Yreceivers = [Ysources(1:22), ones(1,6)*psouth, Ysources(29:end)];
Zreceivers = [ones(1,16)*pbottom, Zsources(17:28), ones(1,4)*pbottom];
% combine all three coordinates
Receivers = [Xreceivers; Yreceivers; Zreceivers];

%% plot sources and receivers
switch debug
    case 'yes'
        figure(f_main);
        % plot source locations
        plot3(Xsources,Ysources,Zsources,'or')
        % plot receiver locations
        hold on
        plot3(Xreceivers,Yreceivers,Zreceivers,'ob')
        % add ellipse
        plot3(Elli1(1,:),Elli1(2,:),Elli1(3,:),'k')
        % fix figure
        axis equal
        xlabel('Easting')
        ylabel('Northing')
        zlabel('Depth')
end

%% compute for series of source-receiver pairs instead
% sources on top %%% what the nubmbers?
SStmp = [5*ones(1,6) 6*ones(1,6) 7*ones(1,6) 8*ones(1,6) 13*ones(1,6)...
    14*ones(1,6) 15*ones(1,6) 16*ones(1,6)];
% receivers on the side %%% numbers?
RRtmp = [[23 24 25 26 27 28], [23 24 25 26 27 28], [23 24 25 26 27 28],...
    [23 24 25 26 27 28], [17 18 19 20 21 22], [17 18 19 20 21 22],...
    [17 18 19 20 21 22], [17 18 19 20 21 22]];
% and vice-versa, combined %%% what?
SS = [SStmp RRtmp];
RR = [RRtmp SStmp];

end
