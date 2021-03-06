function SQR = NoGoSqr(TrialRecord)
% create a plain grey square based on parameters defined in the TrialRecord struct.

if ~isfield(TrialRecord, 'PixelsPerDegree'),
    PixelsPerDegree = 1;
else 
    PixelsPerDegree = TrialRecord.PixelsPerDegree;
end

if ~isfield(TrialRecord, 'sqrsz'),
    sqrsz = round(1 * TrialRecord.PixelsPerDegree);
else 
    sqrsz = round(TrialRecord.sqrsz * TrialRecord.PixelsPerDegree);  %%% ToDo: output the actual size
end

if ~isfield(TrialRecord, 'NGsqrcol'),
    sqrcol = 0.8;
else 
    sqrcol = TrialRecord.NGsqrcol;
end

SQR = repmat(sqrcol, sqrsz, sqrsz);
