function [sta, tcorr] = STA(Iapp, spikes, dt, tminus, tplus)
    if (~exist('tminus'))
        tminus = .075;
    end
    if (~exist('tplus'))
        tplus = .025;
    end
    
    nminus = tminus / dt;               % Sets number of bins below the time of the spike
    nplus = tplus / dt;                 % Sets the number of bins above
    
    tcorr = -nminus*dt:dt:nplus*dt;     % Creates a time window vector
    sta = zeros(size(tcorr));           % Creates spike-triggered average vector of the same size as the time window
    
    spikesum = sum(spikes);             % Calulcates the total number of spikes
    spiketimebins = find(spikes);       % Finds time bins containing spikes and stores them in a vector
    
    ignoredtimebins = 0;                % Will track the number of bins with spikes that are skipped
    
    for j = 1:spikesum                                  % Iterates through the vector containing the times of spikes
        bin = spiketimebins(j);                         % Records the corresponding time
        if (bin > 75 && bin < length(spikes) - 25)      % Tests if the time is greater than 75 and 255 less than the length of spike
            range = (bin - nminus) : (bin + nplus);     % If so, sets a range that begins 75 time bins below and ends 25 bins above the current time 
            sta = sta + Iapp(range);                    % Adds the applied current at each of these time values to the spike-triggered average vector
        else
            ignoredtimebins = ignoredtimebins + 1;      % If not, increases the number of skipped bins by 1
        end
    end
    
    spikesum = spikesum - ignoredtimebins;              % Substracts the total number of spikes by the ignored bins
    sta = sta ./ spikesum;                              % Finds the average
end