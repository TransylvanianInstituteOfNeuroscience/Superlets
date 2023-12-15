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
import superlet 

class SuperletTransformCX:
    """
    Class used to compute the complex Superlet Transform of input data.
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

        self.inputSize          = inputSize
        self.orders             = np.linspace(start=superletOrders[0], stop=superletOrders[1], num=frequencyBins)
        self.convBuffer         = np.zeros(inputSize, dtype=np.complex128)
        self.phaseBuffer        = np.zeros(inputSize, dtype=np.complex128)
        self.absBuffer          = np.zeros(inputSize, dtype=np.float64)
        self.magnitudeBuffer    = np.zeros(inputSize, dtype=np.float64)
        self.superlets          = []

        # create wavelets
        for iFreq in range(frequencyBins):
            centerFreq  = self.frequencies[iFreq]
            nWavelets   = int(np.ceil(self.orders[iFreq]))

            self.superlets.append([])
            for iWave in range(nWavelets):

                # create morlet wavelet
                self.superlets[iFreq].append(superlet.morlet(centerFreq, (iWave + 1) * baseCycles, samplingRate))



    def clear(self):
        """
        Clear the transform.
        """
        # fields
        self.inputSize          = None
        self.superlets          = None
        self.absBuffer          = None
        self.convBuffer         = None
        self.frequencies        = None
        self.orders             = None
        self.magnitudeBuffer    = None
        self.phaseBuffer        = None

    
    def transform(self, inputData):
        """
        Apply the superlet transform on a single data buffer.
        Arguments:
            inputData: A 1xInputSize array containing the signal to be transformed.
        """
        if len(inputData.shape) > 1 or inputData.shape[0] != self.inputSize:
            raise "Invalid input size."

        result = np.zeros((len(self.frequencies), self.inputSize), dtype=np.complex128)

        for iFreq in range(len(self.frequencies)):
            
            # init pooling buffer
            self.magnitudeBuffer.fill(1)
            self.phaseBuffer    .fill(0)

            if len(self.superlets[iFreq]) > 1:
                
                # superlet
                nWavelets   = int(np.floor(self.orders[iFreq]))
                rfactor     = 1.0 / nWavelets

                for iWave in range(nWavelets):
                    self.convBuffer         = fftconvolve(inputData, self.superlets[iFreq][iWave], "same")
                    self.absBuffer          = np.abs(self.convBuffer)
                    self.magnitudeBuffer    *= self.absBuffer
                    self.phaseBuffer        += self.convBuffer / self.absBuffer

                if superlet.fractional(self.orders[iFreq]) != 0 and len(self.superlets[iFreq]) == nWavelets + 1:

                    # apply the fractional wavelet
                    exponent    = self.orders[iFreq] - nWavelets
                    rfactor     = 1 / (nWavelets + exponent)

                    self.convBuffer         = fftconvolve(inputData, self.superlets[iFreq][nWavelets], "same")
                    self.absBuffer          = np.abs(self.convBuffer)
                    self.magnitudeBuffer    *= self.absBuffer ** exponent
                    self.phaseBuffer        += (self.convBuffer / self.absBuffer) * exponent

                # perform geometric mean
                result[iFreq, :] += (self.magnitudeBuffer ** rfactor) * (self.phaseBuffer / self.orders[iFreq])

            else:
                # wavelet transform
                result[iFreq, :] += fftconvolve(inputData, self.superlets[iFreq][0], "same")

        return result


# main superlet function
def superlets(data,
              fs,
              foi,
              c1,
              ord):
    """
    Perform fractional adaptive superlet transform (FASLT) on a single trial.
    Arguments:
        data: a numpy array of data
        fs: the sampling rate in Hz
        foi: list of frequencies of interest
        c1: base number of cycles parameter
        ord: the order (for SLT) or order range (for FASLT), spanned across the frequencies of interest
    Returns: a matrix containing the complex superlet spectrum
    """
    # determine buffer size
    bufferSize = data.shape[len(data.shape) - 1]

    # make order parameter
    if len(ord) == 1:
        ord = (ord, ord)

    # build the superlet analyzer
    faslt = SuperletTransformCX(inputSize        = bufferSize, 
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
