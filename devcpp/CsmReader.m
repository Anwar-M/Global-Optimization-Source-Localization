% Compatible with HDF5 CSM definitions Rev 2.2

h5in_file='Benchmark0.h5';

% READ HDF5 DATA

RevNrMaj=h5readatt(h5in_file,'/MetaData','revisionNumberMajor');
RevNrMin=h5readatt(h5in_file,'/MetaData','revisionNumberMinor');
num_mic=h5readatt(h5in_file,'/MetaData/ArrayAttributes','microphoneCount');
mic_info=h5read(h5in_file,'/MetaData/ArrayAttributes/microphonePositionsM');
coor_ref=h5readatt(h5in_file,'/MetaData/TestAttributes','coordinateReference');
flow_type=h5readatt(h5in_file,'/MetaData/TestAttributes','flowType');
test_description=h5readatt(h5in_file,'/MetaData/TestAttributes','testDescription');
bounds=h5read(h5in_file,'/MetaData/TestAttributes/domainBoundsM');

mach_vector=h5read(h5in_file,'/MeasurementData/machNumber');
rel_hum=h5read(h5in_file,'/MeasurementData/relativeHumidityPct');
csound=h5read(h5in_file,'/MeasurementData/speedOfSoundMPerS');
pres=h5read(h5in_file,'/MeasurementData/staticPressurePa');
temp=h5read(h5in_file,'/MeasurementData/staticTemperatureK');

freqs=h5read(h5in_file,'/CsmData/binCenterFrequenciesHz');
num_freq=h5readatt(h5in_file,'/CsmData/binCenterFrequenciesHz','frequencyBinCount');
csmUnits=h5readatt(h5in_file,'/CsmData','csmUnits');
fftSign=h5readatt(h5in_file,'/CsmData','fftSign');
spectrumType=h5readatt(h5in_file,'/CsmData','spectrumType');

cpimag=h5read(h5in_file,'/CsmData/csmImaginary');
cpreal=h5read(h5in_file,'/CsmData/csmReal');