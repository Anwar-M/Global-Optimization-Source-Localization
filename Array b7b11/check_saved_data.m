clc;
data_path = '.\DATA';
files = dir([data_path '\*.mat']);

tic;
freqs_all = (0:50:25600).';
for I = 1:length(files)
    name = files(I).name;
    for J = 1:length(freqs_all)
        if (str2num(name(5:end-6)))==freqs_all(J)
            freqs_all(J, 2) = 1;
        end
    end
end
fprintf('\tFrequency points processed: %d/%d\n', sum(freqs_all(:,2)), length(freqs_all));

if (sum(freqs_all(:,2))==length(freqs_all))
    fprintf('\tAll frequencies processed!\n');
else
    fprintf('\tProcessed data not complete!\n');
end