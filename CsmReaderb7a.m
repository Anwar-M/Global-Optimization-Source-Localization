h5in_file='ab7aCsmEss.h5';

% READ HDF5 DATA

csound=h5readatt(h5in_file,'/MeasurementData', 'speedOfSoundMPerS');

freqs = h5read(h5in_file, '/CsmData/binCenterFrequenciesHz');
cpimag = h5read(h5in_file, '/CsmData/csmImaginary');
cpreal = h5read(h5in_file, '/CsmData/csmReal');

num_mic=h5readatt(h5in_file,'/MetaData/ArrayAttributes','microphoneCount');
mic_info=h5read(h5in_file,'/MetaData/ArrayAttributes/microphonePositionsM');