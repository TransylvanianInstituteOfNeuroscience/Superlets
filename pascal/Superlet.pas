{
 Project:     Generic
 Module:      Superlet.pas
 Description: A class that computes a superlet with a particular wavelet instantiation and a family of expanded wavelets with the same central frequency
 Author:      Muresan Raul
 History:     18.02.2019 - Created and implemented the first version
}

unit Superlet;

interface

uses Declare,Classes,Math,SysUtils,SingleList,dspComplex,dspTypes,Wavelet,Convolver,ConvolverSuperOptimized,PointerList,DoubleList;

const SL_CONV_SEGMENT_SECONDS = 0.25;
      //SL_CONV_SEGMENT_SECONDS = 1;
      SL_NORMALIZE_WAVELET = TRUE;

type
    //The single superlet manager
    TSuperlet = class
        private
               //Buffers for the real and imaginary part of the spectrum in time
               daMagnitudes:TDoubleList;
               saReal,saImag:TSingleList;
               iNrTimeBins:integer;
               iWaveletType:integer;

               //Wavelets at different scales but same frequency
               paWaveletSet:TPointerList;

               //Convolvers
               paConvolvers:TPointerList;

               //Internal variables
               iNrInputSamples:integer;
               iSampleStepSize:integer;
               sMsStepSize:single;
               sSamplingRate:single;
               sCentralFreqHz:single;
               sNrBaseCycles:single;
               sBandwidthSamples:single;
               bUseExtension:boolean;
               sOrder:single;
               iIntegerOrder:integer;
               sFractionalOrder:single;
               iNrWavelets:integer;
               iNecessaryExtension:integer;
               bBandwidthCompensation:boolean;
               bAdditiveOrder:boolean;

               //Management interface
               function GetSpectralMagnitude(idx:integer):double; inline;
               function GetSpectralPower(idx:integer):double; inline;
               function GetWaveletBaseSize:integer; inline;

               procedure _ComputeDC(psaSourceSignal:PSingleExtArray);

        public
               property SpectralMagnitude[idx:integer]:double read GetSpectralMagnitude;
               property SpectralPower[idx:integer]:double read GetSpectralPower;
               property NrTimeBins:integer read iNrTimeBins;
               property NrInputTimeBins:integer read iNrInputSamples;
               property SamplingRate:single read sSamplingRate;
               property CentralFrequency:single read sCentralFreqHz;
               property StepSizeSamples:integer read iSampleStepSize;
               property StepSizeMs:single read sMsStepSize;
               property WaveletType:integer read iWaveletType;
               property NrBaseCycles:single read sNrBaseCycles;
               property UseExtension:boolean read bUseExtension;
               property WaveletBaseSize:integer read GetWaveletBaseSize;
               property BandwidthSamples:single read sBandwidthSamples;
               property Order:single read sOrder;
               property NecessaryExtension:integer read iNecessaryExtension;
               property AdditiveOrder:boolean read bAdditiveOrder;

               //Construction and destruction interface
               constructor Create(piWaveletType:integer;psCentralFreqHz,psSamplingRateHz,psNrBaseCycles,psOrder:single;piNrInputSamples,piStepSamples:integer;pbAdditiveOrder,pbUseExtension:boolean);
               destructor  Destroy; override;

               //Functional interface
               function  Assign(ppSourceSuperlet:TSuperlet):integer;                                    //Copies another single wavelet spectrum into this one
               function  ComputeSuperlet(psaSourceSignal:PSingleExtArray;pbBandwidthCompensation:boolean):integer; inline; //Computes the super wavelet spectrum of a source signal
               function  GetPower(piTime:integer):double; inline;                                              //Computes the power value associated to a time point
               function  GetPowerLogged(piTime:integer;psLogOffset:single):double; inline;                     //Computes the power value associated to a time point
               function  GetAbs(piTime:integer):double; inline;                                                //Computes the magnitude value associated to a time point
               function  GetAbsLogged(piTime:integer;psLogOffset:single):double; inline;                       //Computes the magnitude value associated to a time point
               function  GetSqrPower(piTime:integer):double; inline;                                           //Computes the squared power value associated to a time point
               function  GetSqrPowerLogged(piTime:integer;psLogOffset:single):double; inline;                  //Computes the squared power value associated to a time point
    end;

implementation

uses Utils,MemoryHandling,Morlet,Windows,FFTW3FloatLibInterface,uTExtendedX87;

//Construction and destruction interface ------------------------------------------------------------------------------------------------------------------------------------
constructor TSuperlet.Create(piWaveletType:integer;psCentralFreqHz,psSamplingRateHz,psNrBaseCycles,psOrder:single;piNrInputSamples,piStepSamples:integer;pbAdditiveOrder,pbUseExtension:boolean);
var iExtension:integer;
    pWavelet:TWavelet;
    pConvolver:TConvolverSuperOptimized;
    iOrder:integer;
    dBandwidth:double;
begin
     sSamplingRate := psSamplingRateHz;
     sCentralFreqHz := psCentralFreqHz;
     sNrBaseCycles := psNrBaseCycles;
     bUseExtension := pbUseExtension;
     sBandwidthSamples := 0;
     sOrder := psOrder;
     iIntegerOrder := trunc(psOrder);
     sFractionalOrder := sOrder - iIntegerOrder;
     if(sFractionalOrder = 0) then iNrWavelets := iIntegerOrder
     else iNrWavelets := iIntegerOrder + 1;
     bAdditiveOrder := pbAdditiveOrder;

     iSampleStepSize := max(piStepSamples,1);
     sMsStepSize := SamplesToMilliseconds(iSampleStepSize,sSamplingRate);

     if(iSampleStepSize = 1) then iNrTimeBins := piNrInputSamples
     else
     begin
          if(piNrInputSamples mod iSampleStepSize = 0) then iNrTimeBins := piNrInputSamples div iSampleStepSize
          else iNrTimeBins := piNrInputSamples div iSampleStepSize + 1;
     end;

     saReal := TSingleList.CreateSized(iNrTimeBins,false);
     saImag := TSingleList.CreateSized(iNrTimeBins,false);
     daMagnitudes := TDoubleList.CreateSized(iNrTimeBins,false);

     paWaveletSet := TPointerList.Create;
     paConvolvers := TPointerList.Create;
     iExtension := 0;
     sBandwidthSamples := 0;
     dBandwidth := 1;

     if(sCentralFreqHz <> 0) then
     begin
          //Create the wavelets
          for iOrder:=1 to iNrWavelets do
          begin
               pWavelet := TWavelet.Create;

               if(not pbAdditiveOrder) then pWavelet.ComputeMorletAuto(sCentralFreqHz,sSamplingRate,sNrBaseCycles * iOrder,SL_NORMALIZE_WAVELET)
               else pWavelet.ComputeMorletAuto(sCentralFreqHz,sSamplingRate,sNrBaseCycles + iOrder - 1,SL_NORMALIZE_WAVELET);

               paWaveletSet.Add(pWavelet);
               dBandwidth := dBandwidth * pWavelet.BandwidthSamples;
               if(pbUseExtension) then
                   if(pWavelet.HalfSize > iExtension) then iExtension := pWavelet.HalfSize;
          end;
          sBandwidthSamples := power(dBandwidth,1/sOrder);

          //Create the convolvers
          for iOrder:=1 to iNrWavelets do
          begin
               pWavelet := TWavelet(paWaveletSet[iOrder-1]);
               pConvolver := TConvolverSuperOptimized.Create(pWavelet.Size,piNrInputSamples,round(SL_CONV_SEGMENT_SECONDS * sSamplingRate),iExtension,
                                                             pWavelet.saRealBuffer,pWavelet.saImagBuffer,saReal.saArray,saImag.saArray);
               paConvolvers.Add(pConvolver);
          end;
     end
     else sBandwidthSamples := iNrInputSamples;

     iNrInputSamples := piNrInputSamples;
     bBandwidthCompensation := false;
     iNecessaryExtension := iExtension;
end;

destructor TSuperlet.Destroy;
var i:integer;
begin
     saReal.Free;
     saImag.Free;
     daMagnitudes.Free;
     for i:=0 to paWaveletSet.Count - 1 do TWavelet(paWaveletSet[i]).Free;
     paWaveletSet.Free;
     for i:=0 to paConvolvers.Count - 1 do TConvolverSuperOptimized(paConvolvers[i]).Free;
     paConvolvers.Free;
end;

//Management interface ------------------------------------------------------------------------------------------------------------------------------------------------------
function TSuperlet.GetSpectralMagnitude(idx:integer):double;
begin
     result := daMagnitudes.daArray[idx];
end;

function TSuperlet.GetSpectralPower(idx:integer):double;
begin
     result := sqr(daMagnitudes.daArray[idx]);
end;

function TSuperlet.GetWaveletBaseSize:integer;
begin
     if(paWaveletSet.Count > 0) then result := TWavelet(paWaveletSet[0]).Size
     else result := iNrInputSamples;
end;

procedure TSuperlet._ComputeDC(psaSourceSignal:PSingleExtArray);
var i:integer;
    dSum:double;
    sAvg:single;
begin
     sAvg := 0;
     dSum := 0;
     for i:=0 to iNrInputSamples - 1 do dSum := dSum + psaSourceSignal[i];
     if(iNrInputSamples <> 0) then sAvg := dSum / iNrInputSamples;

     saImag.InitializeZero;
     saReal.Initialize(sAvg);
     daMagnitudes.Initialize(sAvg);
end;

//Functional interface ------------------------------------------------------------------------------------------------------------------------------------------------------
//Copies another convolver list into this one
function TSuperlet.Assign(ppSourceSuperlet:TSuperlet):integer;
var i:integer;
begin
     result := -1;
     if((ppSourceSuperlet = nil) or (ppSourceSuperlet.sOrder <> sOrder)) then exit;

     saReal.Assign(ppSourceSuperlet.saReal);
     saImag.Assign(ppSourceSuperlet.saImag);
     daMagnitudes.Assign(ppSourceSuperlet.daMagnitudes);

     for i:=0 to paWaveletSet.Count - 1 do TWavelet(paWaveletSet[i]).Assign(TWavelet(ppSourceSuperlet.paWaveletSet[i]));
     for i:=0 to paConvolvers.Count - 1 do TConvolverSuperOptimized(paConvolvers[i]).Assign(TConvolverSuperOptimized(ppSourceSuperlet.paConvolvers[i]));

     iNrTimeBins            := ppSourceSuperlet.iNrTimeBins;
     iSampleStepSize        := ppSourceSuperlet.iSampleStepSize;
     sMsStepSize            := ppSourceSuperlet.sMsStepSize;
     sSamplingRate          := ppSourceSuperlet.sSamplingRate;
     sCentralFreqHz         := ppSourceSuperlet.sCentralFreqHz;
     sBandwidthSamples      := ppSourceSuperlet.sBandwidthSamples;
     bUseExtension          := ppSourceSuperlet.bUseExtension;
     sNrBaseCycles          := ppSourceSuperlet.sNrBaseCycles;
     sOrder                 := ppSourceSuperlet.sOrder;
     iIntegerOrder          := ppSourceSuperlet.iIntegerOrder;
     sFractionalOrder       := ppSourceSuperlet.sFractionalOrder;
     iNrWavelets            := ppSourceSuperlet.iNrWavelets;
     iNecessaryExtension    := ppSourceSuperlet.iNecessaryExtension;
     bBandwidthCompensation := ppSourceSuperlet.bBandwidthCompensation;
     bAdditiveOrder         := ppSourceSuperlet.bAdditiveOrder;

     result := 0;
end;

//Computes the superlet of a source signal
function TSuperlet.ComputeSuperlet(psaSourceSignal:PSingleExtArray;pbBandwidthCompensation:boolean):integer;
var iResult:integer;
    sBandwidthCompensationFactor:single;
    iWaveletIdx,iTime:integer;
    dScaling:double;
    bForceClassicConvolution:boolean;
begin
     result := -1;

     if(psaSourceSignal = nil) then exit;
     bBandwidthCompensation := pbBandwidthCompensation;


     if(sCentralFreqHz <> 0) then
     begin
          //Initialize the magnitude buffer with 1
          daMagnitudes.InitializeOne;
          dScaling := sqrt(2);

          for iWaveletIdx:=1 to iNrWavelets do
          begin
               //Force classic convolution for the shortest wavelet to avoid FFT-based convolution issues with small kernels
               if(iWaveletIdx = 1) then bForceClassicConvolution := true
               else bForceClassicConvolution := false;

               //Compute the convolution
               if(iSampleStepSize = 1) then iResult := TConvolverSuperOptimized(paConvolvers[iWaveletIdx-1]).ConvolveOA(psaSourceSignal,bForceClassicConvolution)
               else iResult := TConvolverSuperOptimized(paConvolvers[iWaveletIdx-1]).ConvolveOA_Stepped(psaSourceSignal,iSampleStepSize,bForceClassicConvolution);
               if(iResult <> 0) then exit(-2);

               //Compensate bandwidth
               if(bBandwidthCompensation) then
               begin
                    sBandwidthCompensationFactor := Math.Max(sqrt(sCentralFreqHz / sNrBaseCycles * WP_MORLET_SD_SPREAD),1.0);
                    saReal.DivScalar(sBandwidthCompensationFactor);
                    saImag.DivScalar(sBandwidthCompensationFactor);
               end;

               //Compute the geometric product, taking care of the fractional wavelet case
               if((iWaveletIdx < iNrWavelets) or (sFractionalOrder = 0)) then for iTime:=0 to iNrTimeBins - 1 do daMagnitudes.daArray[iTime] := daMagnitudes.daArray[iTime] * sqrt(sqr(saReal.saArray[iTime])+sqr(saImag.saArray[iTime]))
               else for iTime:=0 to iNrTimeBins - 1 do daMagnitudes.daArray[iTime] := daMagnitudes.daArray[iTime] * power(sqrt(sqr(saReal.saArray[iTime])+sqr(saImag.saArray[iTime])),sFractionalOrder);
          end;

          //Compute the geometric mean
          for iTime:=0 to iNrTimeBins - 1 do daMagnitudes.daArray[iTime] := dScaling * power(daMagnitudes.daArray[iTime],1/sOrder);
     end
     else _ComputeDC(psaSourceSignal);

     result := 0;
end;

//Computes the power value associated to a time point
function TSuperlet.GetPower(piTime:integer):double;
begin
     result := sqr(daMagnitudes[piTime]);
end;

//Computes the power value associated to a time point
function TSuperlet.GetPowerLogged(piTime:integer;psLogOffset:single):double;
begin
     result := log10(psLogOffset + sqr(daMagnitudes[piTime]));
end;

//Computes the magnitude value associated to a time point
function TSuperlet.GetAbs(piTime:integer):double;
begin
     result := daMagnitudes[piTime];
end;

//Computes the magnitude value associated to a time point
function TSuperlet.GetAbsLogged(piTime:integer;psLogOffset:single):double;
begin
     result := log10(psLogOffset + daMagnitudes[piTime]);
end;

//Computes the squared power value associated to a time point
function TSuperlet.GetSqrPower(piTime:integer):double;
begin
     result := sqr(sqr(daMagnitudes[piTime]));
end;

//Computes the squared power value associated to a time point
function TSuperlet.GetSqrPowerLogged(piTime:integer;psLogOffset:single):double;
begin
     result := log10(psLogOffset + sqr(sqr(daMagnitudes[piTime])));
end;


end.
