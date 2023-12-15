# Time-frequency analysis with superlets
# Based on 'Time-frequency super-resolution with superlets'
# by Moca et al., 2021 Nature Communications
#
# Implementation by Harald Bârzan and Richard Eugen Ardelean

#
# Note: for runs on multiple batches of data, the class SuperletTransform can be instantiated just once
# this saves time and memory allocation for the wavelets and buffers
#


import numpy as np
from scipy.signal import fftconvolve

# spread, in units of standard deviation, of the Gaussian window of the Morlet wavelet
MORLET_SD_SPREAD = 6

# the length, in units of standard deviation, of the actual support window of the Morlet
MORLET_SD_FACTOR = 2.5



def computeWaveletSize(fc, nc, fs):
    """
    Compute the size in samples of a morlet wavelet.
    Arguments:
        fc - center frequency in Hz
        nc - number of cycles
        fs - sampling rate in Hz
    """
    sd = (nc / 2) * (1 / np.abs(fc)) / MORLET_SD_FACTOR
    return int(2 * np.floor(np.round(sd * fs * MORLET_SD_SPREAD) / 2) + 1)


def gausswin(size, alpha):
    """
    Create a Gaussian window.
    """
    halfSize    = int(np.floor(size / 2))
    idiv        = alpha / halfSize

    t = (np.array(range(size), dtype=np.float64) - halfSize) * idiv
    window = np.exp(-(t * t) * 0.5)
    
    return window

    

def morlet(fc, nc, fs):
    """
    Create an analytic Morlet wavelet.
    Arguments:
        fc - center frequency in Hz
        nc - number of cycles
        fs - sampling rate in Hz
    """
    size    = computeWaveletSize(fc, nc, fs)
    half    = int(np.floor(size / 2))
    gauss   = gausswin(size, MORLET_SD_SPREAD / 2)
    igsum   = 1 / gauss.sum()
    ifs     = 1 / fs

    t = (np.array(range(size), dtype=np.float64) - half) * ifs
    wavelet = gauss * np.exp(2 * np.pi * fc * t * 1j) * igsum

    return wavelet

def fractional(x):
    """
    Get the fractional part of the scalar value x.
    """
    return x - int(x)


class SuperletTransform:
    """
    Class used to compute the Superlet Transform of input data.
    """

    def __init__(   self,
                    inputSize,
                    samplingRate,
                    frequencyRange,
                    frequencyBins,
                    baseCycles,
                    superletOrders,
                    frequencies = None):
        """
        Initialize the superlet transform. 
        Arguments:
            inputSize: size of the input in samples
            samplingRate: the sampling rate of the input signal in Hz
            frequencyRange: tuplet of ascending frequency points, in Hz
            frequencyBins: number of frequency bins to sample in the interval frequencyRange
            baseCycles: number of cycles of the smallest wavelet (c1 in the paper)
            superletOrders: a tuple containing the range of superlet orders, linearly distributed along frequencyRange
            frequencies: specific list of frequencies - can be provided in stead of frequencyRange (it is ignored in this case)
        """
        # clear to reinit
        self.clear()

        # initialize containers
        if frequencies is not None:
            frequencyBins = len(frequencies)
            self.frequencies = frequencies
        else:
            self.frequencies = np.linspace(start=frequencyRange[0], stop=frequencyRange[1], num=frequencyBins)

        self.inputSize      = inputSize
        self.orders         = np.linspace(start=superletOrders[0], stop=superletOrders[1], num=frequencyBins)
        self.convBuffer     = np.zeros(inputSize, dtype=np.complex128)
        self.poolBuffer     = np.zeros(inputSize, dtype=np.float64)
        self.superlets      = []

        # create wavelets
        for iFreq in range(frequencyBins):
            centerFreq  = self.frequencies[iFreq]
            nWavelets   = int(np.ceil(self.orders[iFreq]))

            self.superlets.append([])
            for iWave in range(nWavelets):

                # create morlet wavelet
                self.superlets[iFreq].append(morlet(centerFreq, (iWave + 1) * baseCycles, samplingRate))


    def __del__(self):
        """
        Destructor.
        """
        self.clear()


    def clear(self):
        """
        Clear the transform.
        """
        # fields
        self.inputSize   = None
        self.superlets   = None
        self.poolBuffer  = None
        self.convBuffer  = None
        self.frequencies = None
        self.orders      = None


    
    def transform(self, inputData):
        """
        Apply the transform to a buffer or list of buffers.
        Arguments:
            inputData - an NDarray of input data
        """

        # compute number of arrays to transform
        if len(inputData.shape) == 1:
            if inputData.shape[0] != self.inputSize:
                raise "Input data must meet the defined input size for this transform."
            
            result = np.zeros((self.inputSize, len(self.frequencies)), dtype=np.float64)
            self.transformOne(inputData, result)
            return result

        else:
            n       = int(np.sum(inputData.shape[0:len(inputData.shape) - 1]))
            insize  = int(inputData.shape[len(inputData.shape) - 1])

            if insize != self.inputSize:
                raise "Input data must meet the defined input size for this transform."
            
            # reshape to data list
            datalist = np.reshape(inputData, (n, insize), 'C')
            result = np.zeros((len(self.frequencies), self.inputSize), dtype=np.float64)

            for i in range(0, n):
                self.transformOne(datalist[i, :], result)

            return result / n


    def transformOne(self, inputData, accumulator):
        """
        Apply the superlet transform on a single data buffer.
        Arguments:
            inputData: A 1xInputSize array containing the signal to be transformed.
            accumulator: a spectrum to accumulate the resulting superlet transform
        """
        accumulator.resize((len(self.frequencies), self.inputSize))

        for iFreq in range(len(self.frequencies)):
            
            # init pooling buffer
            self.poolBuffer.fill(1)

            if len(self.superlets[iFreq]) > 1:
                
                # superlet
                nWavelets   = int(np.floor(self.orders[iFreq]))
                rfactor     = 1.0 / nWavelets

                for iWave in range(nWavelets):
                    self.convBuffer = fftconvolve(inputData, self.superlets[iFreq][iWave], "same")
                    self.poolBuffer *= 2 * np.abs(self.convBuffer) ** 2

                if fractional(self.orders[iFreq]) != 0 and len(self.superlets[iFreq]) == nWavelets + 1:

                    # apply the fractional wavelet
                    exponent = self.orders[iFreq] - nWavelets
                    rfactor = 1 / (nWavelets + exponent)

                    self.convBuffer = fftconvolve(inputData, self.superlets[iFreq][nWavelets], "same")
                    self.poolBuffer *= (2 * np.abs(self.convBuffer) ** 2) ** exponent

                # perform geometric mean
                accumulator[iFreq, :] += self.poolBuffer ** rfactor


            else:
                # wavelet transform
                accumulator[iFreq, :] += (2 * np.abs(fftconvolve(inputData, self.superlets[iFreq][0], "same")) ** 2).astype(np.float64)


# main superlet function
def superlets(data,
              fs,
              foi,
              c1,
              ord):
    """
    Perform fractional adaptive superlet transform (FASLT) on a list of trials. 
    Arguments:
        data: a numpy array of data. The rightmost dimension of the data is the trial size. The result will be the average over all the spectra.
        fs: the sampling rate in Hz
        foi: list of frequencies of interest
        c1: base number of cycles parameter
        ord: the order (for SLT) or order range (for FASLT), spanned across the frequencies of interest
    Returns: a matrix containing the average superlet spectrum
    """
    # determine buffer size
    bufferSize = data.shape[len(data.shape) - 1]

    # make order parameter
    if len(ord) == 1:
        ord = (ord, ord)

    # build the superlet analyzer
    faslt = SuperletTransform(  inputSize        = bufferSize, 
                                frequencyRange   = None, 
                                frequencyBins    = None, 
                                samplingRate     = fs, 
                                frequencies      = foi, 
                                baseCycles       = c1, 
                                superletOrders   = ord)
        
    # apply transform
    result = faslt.transform(data)
    faslt.clear()

    return result
