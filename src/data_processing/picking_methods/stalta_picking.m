% Dong Liu adds this part on 12/09/2019
% this is to pick the arrival diffraction for several sequences between two
% manually picking points. By using this function a few times, we can
% have all the interpolated lines and its best-fit correlated arrival time
% Here the picking strategy is the short-term vs long-term average
function [tlist]=stalta_picking(dataSubs,windowsize,dt,xse,yse,endnoise)
% windowsize=half window size
% for this specif function it's better to have relatively large window size
% the correlation is done with the starting sequence
xs=xse(1);
xe=xse(2);
ys=yse(1);
ye=yse(2);

x=xs:xe;
ymidlist=floor((ys+(x-xs)/(xe-xs)*(ye-ys))+0.5);
ylowlist=ymidlist-windowsize;
yuplist=ymidlist+windowsize;
tneighbor=zeros(1,size(x,2));
for it=1:xe-xs+1
    tindex=arrival_time(dataSubs(xs+it-1,1:end));
    % need to replace the cross-correlation with stalta picking method    
    if tindex==0
        tneighbor(1,it)=ymidlist(1,it)*dt*10^6;
    else
        tneighbor(1,it)=(endnoise-1+tindex)*dt*10^6; % arrival time (unit in us) this is because the input signal comes from the removed noise signals
    end
    
% % one can use this part to control the range of the index    
%      if tneighbor(1,it)<ylowlist(1,it)*dt*10^6
%          tneighbor(1,it)=ymidlist(1,it)*dt*10^6;
%      elseif tneighbor(1,it)>yuplist(1,it)*dt*10^6
%          tneighbor(1,it)=ymidlist(1,it)*dt*10^6;
%      end
end 
tlist=tneighbor;
end


% for one signal, get the most suitable arrival time
function tindex=arrival_time(data1seq)
abs_v = abs(data1seq);% it should be a 1xN array
l_v=size(data1seq,2);% length of the datapoint
trig_array = zeros(1,2);% initialization
ntrig=0;
lta_calc_flag = 0; % to test if the full calculation is needed or not?

% set all the controlling parameters
l_sta = round(4/160*8000);     % STA window length 4 us
l_lta = round(20/160*8000);     % LTA window length 20 us
th_on = 1.2;        % Trigger on when sta_to_lta exceeds this theshold
th_off = 1.;     % Trigger off when sta_to_lta drops below threshold
min_dur = round(0/160*8000);   % Any triggers shorter than min_dur are discarded

ratio_trig = zeros(l_v-l_lta,1); % record the stalta_ratio for each i value

% i is the primary reference point (right end of STA/LTA window)
i = l_lta+1;
while i <= l_v % START STA_LTA MAIN LOOP
    if (lta_calc_flag == 0)
      lta_sum = 0;
      sta_sum = 0;
      for j = i-l_lta:i-1              % Loop to compute LTA & STA
         lta_sum = lta_sum + abs_v(j); % Sum LTA window
         if (i - j) <= l_sta           % Sum STA window (right side of LTA)
            sta_sum = sta_sum + abs_v(j);
         end
      end
      lta_calc_flag = 1; % after the first calculation, the following is just a moving window, no need to recalculate everything again
    else
      lta_sum = lta_sum - abs_v(i-l_lta-1) + abs_v(i-1);
      sta_sum = sta_sum - abs_v(i-l_sta-1) + abs_v(i-1);
    end
    lta = lta_sum/l_lta;
    sta = sta_sum/l_sta;
    sta_to_lta = sta/lta;
    
   if (sta_to_lta > th_on)
      j = i;   % Set secondary reference point = primary
      g = 0;   % l_lta growth, only used if LTA growing
      while (sta_to_lta > th_off)
         j = j+1;
         if j < l_v
            sta_sum = sta_sum - abs_v(j-l_sta-1) + abs_v(j-1);
            
%             switch lta_mode
%                case 'frozen'
%                   % LTA is good just the way it is
%                case 'continuous'
            % Add new data point, remove oldest data point
            % we assume here the windowsize of LTA moves together with the
            % STA  to the right
             lta_sum = lta_sum - abs_v(j-l_lta-1) + abs_v(j-1);
%                case 'grow'
%                   % Add new data point, increase
%                   lta_sum = lta_sum + abs_v(j-1);
%                   l_lta = l_lta + 1;
%                   g = g+1;
%             end
            sta = sta_sum/l_sta;
            lta = lta_sum/l_lta;
            sta_to_lta = sta/lta;
         else
            sta_to_lta = 0; % Force trigger off (end of data)
         end
      end
      duration = (j-i); % span from trigger on to trigger off
      l_lta = l_lta-g;
      if duration >= min_dur % If duration < min_dur then skip it
         trig_t = i-l_sta;  % Beginning of STA window during trigger on
         end_t  = j;        % End of STA window during trigger off
         ntrig = ntrig + 1; % Event counter
         trig_array(ntrig,:) = [trig_t, end_t];
      end
      lta_calc_flag = 0;
   end
   ratio_trig(i-l_lta,1)=sta_to_lta;
   i = i + 1;
end

if (trig_array(1,1)==0)&&(trig_array(1,2)==0)
    disp('No events detected')
    tindex = 0;
    return
end

tindex = trig_array(1,1);
disp(tindex);
disp(["t index is ", num2str(tindex), num2str(sta_to_lta)]);
figure
plot((1:size(ratio_trig,1))'/8000*160+l_lta*160/8000,ratio_trig)
xlabel('Travel time (\mu s)')
ylabel('STA/LTA')

end
