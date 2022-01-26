{
 Project:     Generic
 Module:      Superlet.pas
 Description: A class that computes a the Morlet wavelet
 Author:      Muresan Raul
 History:     18.02.2019 - Created and implemented the first version
}

unit Morlet;

interface

uses dspComplex,dspTypes,Math,uTExtendedX87;

         //Computes the general Morlet form (simplified, when B >= 5)
         function MorletGeneral(t,w0:TdspFloat):TDspComplex;
         //Computes the parametrized Morlet
         function MorletBCf(t,pdBandwidth,pdCenterFrenquency:TdspFloat):TDspComplex;
         //Computes the parametrized Morlet with zero mean (admissible - generalized Morlet)
         function MorletBCfGeneralized(t,pdBandwidth,pdCenterFrenquency:TdspFloat):TDspComplex;
         //Computes the Gaussian associated to the parametrized Morlet
         function MorletGauss(t,pdBandwidth:TdspFloat):TDspFloat;

implementation

//Computes the general Morlet form (simplified, when B >= 5)
function MorletGeneral(t,w0:TdspFloat):TDspComplex;
var OneOverPIConst,exp1,exp2,cSigma:TDspFloat;
    cmpxExp1:TDspComplex;
begin
     cSigma := sqrt(1+exp(-sqr(w0))-2*exp(-3/4*sqr(w0)));
     OneOverPIConst := power(PI,-1/4) / cSigma;
     exp1 := -exp(-sqr(w0)/2);
     exp2 := exp(-sqr(t)/2);
     cmpxExp1 := cmpxInit(exp1,0);
     result := cmpxMul(cmpxSub(cmpxExpj(w0*t),cmpxExp1),OneOverPIConst * exp2);
end;

//Computes the parametrized Morlet
function MorletBCf(t,pdBandwidth,pdCenterFrenquency:TdspFloat):TDspComplex;
var cNorm,exp1:TDspFloat;
    cmpxExp1:TDspComplex;
begin
     cNorm := 1/(pdBandwidth * sqrt(2*pi)); // / (182.4377 * sqrt(30));
     exp1 := cNorm * exp(-sqr(t)/(2*sqr(pdBandwidth)));
     cmpxExp1 := cmpxInit(exp1,0);
     result := cmpxMul(cmpxExpj(2*PI*pdCenterFrenquency*t),cmpxExp1);

end;

//Computes the parametrized Morlet with zero mean (admissible - generalized Morlet)
function MorletBCfGeneralized(t,pdBandwidth,pdCenterFrenquency:TdspFloat):TDspComplex;
var cNorm,exp1,cMean:TDspFloat;
    cmpxExp1:TDspComplex;
    cmpxMean:TDspComplex;
begin
     cMean := exp(-sqr(pdBandwidth)/2*sqr(2*pi*pdCenterFrenquency));
     cNorm := 1/(pdBandwidth * sqrt(2*pi));
     exp1 := cNorm * exp(-sqr(t)/(2*sqr(pdBandwidth)));
     cmpxExp1 := cmpxInit(exp1,0);
     cmpxMean := cmpxInit(cMean,0);

     result := cmpxSub(cmpxMul(cmpxExpj(2*PI*pdCenterFrenquency*t),cmpxExp1),cmpxMean);
end;

//Computes the Gaussian associated to the parametrized Morlet
function MorletGauss(t,pdBandwidth:TdspFloat):TDspFloat;
var cNorm,exp1:TDspFloat;
begin
     cNorm := 1/(pdBandwidth * sqrt(2*pi));
     exp1 := cNorm * exp(-sqr(t)/(2*sqr(pdBandwidth)));
     result := exp1;
end;

end.
