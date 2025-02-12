function [cleanPhase filteredTrace Magnitude]= get_phase(trace, fs, Band)
    % calculateThetaPhase - Calculate the theta phase of a signal
    %
    % Syntax: phase = calculateThetaPhase(trace, fs, thetaBand)
    %
      % Ensure trace is a column vector
    trace = trace(:);


    % Design a bandpass filter for the theta band
    [b, a] = butter(2, Band / (fs / 2), 'bandpass');

    % Apply the filter to the valid portion of the trace
    filteredTrace = filtfilt(b, a, trace);

    % Calculate the analytic signal using the Hilbert transform
    analyticSignal = hilbert(filteredTrace);

    % Extract the instantaneous phase
    cleanPhase = angle(analyticSignal);
    filteredTrace=filteredTrace';
    Magnitude =abs(analyticSignal)';

end