function [T_r_char, T_r_event] = Matfun_Tr_DMCA_event(rain, flow, T_r_max, pth)
%% Matlab function to calculate event-based catchment response time 
% Step i)-vii) follows the calculations steps described in [].
% Step i) includes the DMCA-based method introduced by Giani et al. (2021).
% The original code for the DMCA-based method is available: 
% https://github.com/giuliagiani/Tr_DMCA

%% Step i)
% after Giani et al. (2021) DOI: 10.1029/2020WR028201
% DMCA-based methodology to estimate catchment response time:
% this function finds the catchment response time given a rainfall 
% and a streamflow timeseries aligned in time and with the same
% length and time step.

% INPUT
% max_window: Maximum window tested. Set it sensibly according to the resolution of your data (e.g. hourly data, max_window=300 means that time of concentration can be maximum 300hours/2 = 150hours =~ 6days)
% rain: your rainfall timeserie as a row vector 1 x n
% flow: your streamflow timeseries as a row vector 1 x n

% ROUTINE
rain_int = cumsum(rain, 'omitnan');     %cumulating rainfall timeseries (Eq.1)
flow_int = cumsum(flow, 'omitnan');     %cumulating streamflow timeseries (Eq.2)
T = length(rain);                       %length of the timeseries

rain_mean = zeros(floor(T_r_max/2),length(rain)); flow_mean = rain_mean;
flutt_rain = rain_mean; flutt_flow = rain_mean;

F_rain = zeros(1,floor(T_r_max/2)); F_flow = F_rain; F_rain_flow = F_rain; rho = F_rain;
for window = 3:2:T_r_max
    rain_mean((window-1)/2,:) = movmean(rain_int, window); %moving average on the integrated rainfall timeseries (Eq.5)
    flow_mean((window-1)/2,:) = movmean(flow_int, window); %moving average on the integrated streamflow timeseries (Eq.6)
    flutt_rain((window-1)/2,:) = rain_int-rain_mean((window-1)/2,:);
    F_rain((window-1)/2) = (1/(T-window+1))*nansum((flutt_rain((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))).^2); %Squared rainfall fluctuations (Eq.3)
    flutt_flow((window-1)/2,:) = flow_int-flow_mean((window-1)/2,:);
    F_flow((window-1)/2) = (1/(T-window+1))*nansum((flutt_flow((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))).^2); %Squared streamflow fluctuations (Eq.4)
    F_rain_flow((window-1)/2) = (1/(T-window+1))*nansum(flutt_rain((window-1)/2,window-0.5*(window-1):T-0.5*(window-1)).*flutt_flow((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))); %Bivariate rainfall-streamflow fluctuations (Eq.7)
    rho((window-1)/2) = F_rain_flow((window-1)/2)/(sqrt(F_rain((window-1)/2))*sqrt(F_flow((window-1)/2))); %DMCA-based correlation coefficent (Eq.8)
end

% OUTPUT
position_minimum = find(rho == nanmin(rho));
T_r_char = position_minimum;

% DCMA-based event selection method
% INPUT
% thr: threshold value for event selection, acceptable range: 0-0.5,
% suggested range: 0.01-0.05

%% step ii)
flutt_rain_norm = flutt_rain(T_r_char,:)./max(abs(flutt_rain(T_r_char,:)));
flutt_flow_norm = flutt_flow(T_r_char,:)./max(abs(flutt_flow(T_r_char,:)));

%% step iii)
rr_id = find((flutt_rain_norm > quantile(flutt_rain_norm,1-pth))...
    & (quantile(flutt_flow_norm,pth) > flutt_flow_norm));

%% step iv)
rr_did = rr_id(2:end) - rr_id(1:end-1);
rr_sid = rr_id(find(rr_did ~= 1)+1);
rr_sid(T_r_max >= rr_sid) = []; rr_sid(rr_sid >= (numel(flow) - T_r_max)) = [];

%% step v)
T_r_event = zeros(numel(rr_sid),1); rhomin = T_r_event;
for i = 1:numel(rr_sid)
    sid = rr_sid(i)-floor(T_r_max/2);
    eid = rr_sid(i)+T_r_max;
    % repeat the DMCA routine for every event
        % ROUTINE
        rain_int = cumsum(rain(sid:eid), 'omitnan');     %cumulating rainfall timeseries (Eq.1)
        flow_int = cumsum(flow(sid:eid), 'omitnan');     %cumulating streamflow timeseries (Eq.2)
        T = length(rain(sid:eid));                       %length of the timeseries

        rain_mean = zeros(floor(T_r_max/2),length(rain(sid:eid))); flow_mean = rain_mean;
        flutt_rain = rain_mean; flutt_flow = rain_mean;

        F_rain = zeros(1,floor(T_r_max/2)); F_flow = F_rain; F_rain_flow = F_rain; rho = F_rain;
        for window = 3:2:T_r_max
            rain_mean((window-1)/2,:) = movmean(rain_int, window); %moving average on the integrated rainfall timeseries (Eq.5)
            flow_mean((window-1)/2,:) = movmean(flow_int, window); %moving average on the integrated streamflow timeseries (Eq.6)
            flutt_rain((window-1)/2,:) = rain_int-rain_mean((window-1)/2,:);
            F_rain((window-1)/2) = (1/(T-window+1))*nansum((flutt_rain((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))).^2); %Squared rainfall fluctuations (Eq.3)
            flutt_flow((window-1)/2,:) = flow_int-flow_mean((window-1)/2,:);
            F_flow((window-1)/2) = (1/(T-window+1))*nansum((flutt_flow((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))).^2); %Squared streamflow fluctuations (Eq.4)
            F_rain_flow((window-1)/2) = (1/(T-window+1))*nansum(flutt_rain((window-1)/2,window-0.5*(window-1):T-0.5*(window-1)).*flutt_flow((window-1)/2,window-0.5*(window-1):T-0.5*(window-1))); %Bivariate rainfall-streamflow fluctuations (Eq.7)
            rhoevent((window-1)/2) = F_rain_flow((window-1)/2)/(sqrt(F_rain((window-1)/2))*sqrt(F_flow((window-1)/2))); %DMCA-based correlation coefficent (Eq.8)
        end

        % OUTPUT
        position_minimum = find(rhoevent == nanmin(rhoevent));
        trevent = position_minimum;
    
    rhomin(i,1) = min(rhoevent);
    if isempty(trevent)
        T_r_event(i,1) = NaN;
    else
        T_r_event(i,1) = trevent(1);
    end
end

%% step vi)
T_r_event(rhomin >= 0) = [];
T_r_event(T_r_event == ceil(T_r_max/2)-1) = [];

%% step vii)
T_r_event = rmoutliers(T_r_event);
