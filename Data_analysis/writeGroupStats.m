function writeGroupStats(wormdata, index, settings)
%writeGroupStats will iterate through combined wormdata files and aggregate
%data into an excel file
%   Detailed explanation goes here

[pdir, cdir, ~] = fileparts(settings.workingDir);
details.outputname = fullfile(pdir, cdir, 'groupStats.xlsx');
details.genotype = wormdata(1).genotype;
details.column = char(index+'A'-1);

if isfield(wormdata, 'velocity')
    velocity = horzcat(wormdata(:).velocity);
    writeData(velocity, 'Interval (sec)', details)
end

pkTraces = mean(horzcat(wormdata(:).peakTraces), 2, 'omitmissing');

writeData(wormdata(1).intervalVector, 'Interval (sec)', details)
writeData(wormdata(1).amplitudeVector, 'Peak Amplitude (a.u)', details)
writeData(pkTraces, 'Mean Peak Profile (sec)', details)
writeData(wormdata(1).riseVector, 'Rise Time (sec)', details)
writeData(wormdata(1).tauVector, 'Decay Constant (sec)', details)



    function writeData(data, sheetName, details)
        writematrix(details.genotype, details.outputname, Sheet=sheetName, Range=[details.column '1'])
        writematrix(data, details.outputname, Sheet=sheetName, Range=[details.column '2'])
    end

end