path_to_save = 'O:\Array Benchmark\EnergyFunction\DiffEvo\';
load('TestData.mat');

fileID = fopen([path_to_save 'csmreal.bin'],'w');
fwrite(fileID, real(CSM), 'double');
fclose(fileID);
fileID = fopen([path_to_save 'csmimag.bin'],'w');
fwrite(fileID, imag(CSM), 'double');
fclose(fileID);
fileID = fopen([path_to_save 'micpos.bin'],'w');
fwrite(fileID, mic_info(:), 'double');
fclose(fileID);