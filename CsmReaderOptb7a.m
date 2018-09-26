h5in_file='ab7aCsmOpt.h5';

% READ HDF5 DATA

freqs = h5read(h5in_file, '/RefData/binCenterFrequenciesHz');

source0nb = h5read(h5in_file, '/RefData/Source0/narrowBandReference');
source1nb = h5read(h5in_file, '/RefData/Source1/narrowBandReference');
source2nb = h5read(h5in_file, '/RefData/Source2/narrowBandReference');
source3nb = h5read(h5in_file, '/RefData/Source3/narrowBandReference');

source03rd = h5read(h5in_file, '/RefData/Source0/thirdOctBandReference');
source13rd = h5read(h5in_file, '/RefData/Source1/thirdOctBandReference');
source23rd = h5read(h5in_file, '/RefData/Source2/thirdOctBandReference');
source33rd = h5read(h5in_file, '/RefData/Source3/thirdOctBandReference');

source0nb(:,2) = 10*log10(source0nb(:,1)/2/2e-5/2e-5);
source1nb(:,2) = 10*log10(source1nb(:,1)/2/2e-5/2e-5);
source2nb(:,2) = 10*log10(source2nb(:,1)/2/2e-5/2e-5);
source3nb(:,2) = 10*log10(source3nb(:,1)/2/2e-5/2e-5);

source03rd(:,2) = 10*log10(source03rd(:,1)/2/2e-5/2e-5);
source13rd(:,2) = 10*log10(source13rd(:,1)/2/2e-5/2e-5);
source23rd(:,2) = 10*log10(source23rd(:,1)/2/2e-5/2e-5);
source33rd(:,2) = 10*log10(source33rd(:,1)/2/2e-5/2e-5);