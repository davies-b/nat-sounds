# nat-sounds
'nat-sound' is a repository to study decompositions of sounds with a biomimetic approach based on a cochlea-like metamaterial.

Please cite the following reference when using this source:

Ammari H & Davies B. 2020 *A biomimetic basis for auditory processing and the perception of natural sounds* arXiv:2005.12794 (https://arxiv.org/abs/2005.12794)

# Installation

Add the repository 'nat-sound' to the MATLAB file path.

The routine 'RUN_perceptionvec' is used to generate the vector of natural sound coefficients for a given sound sample.

# Running examples

Some examples of audio files are included in the folder 'nat-sounds/samples'. The desired audio file can be selected with the 'file' variable in line.

There is also the option to bypass the use of 'audiofilter' and load pre-filtered data for a recording of a trumpet, stored in the file 'trumpet.mat'.

# Updates

New version (March 2021): code updated to use the 'fconv' routine to compute convolutions efficiently. This achieves a speed up by using the fft/ifft. This is based on a code by Stephen McGovern:

Stephen McGovern (2021). *Fast Convolution* (https://www.mathworks.com/matlabcentral/fileexchange/5110-fast-convolution), MATLAB Central File Exchange. Retrieved March 16, 2021.
