unit Declare;

interface

         type THandle = integer;

         type PPointerExtArray = ^TPointerExtArray;
              TPointerExtArray =  array[0..0] of Pointer;

         type PByteExtArray = ^TByteExtArray;
              TByteExtArray =  array[0..9999] of byte;

         type PIntegerExtArray = ^TIntegerExtArray;
              TIntegerExtArray =  array[0..0] of integer;

         type PInteger16ExtArray = ^TInteger16ExtArray;
              TInteger16ExtArray =  array[0..0] of SmallInt;

         type PIntegerExtArrayArray = ^TIntegerExtArrayArray;
              TIntegerExtArrayArray =  array[0..0] of PIntegerExtArray;

         type PIntegerExtArrayArrayArray = ^TIntegerExtArrayArrayArray;
              TIntegerExtArrayArrayArray =  array[0..0] of PIntegerExtArrayArray;

         type PLongWordExtArray = ^TLongWordExtArray;
              TLongWordExtArray =  array[0..0] of longword;

         type PSingleExtArray  = ^TSingleExtArray;
              TSingleExtArray  =  array[0..0] of single;

         type PSingleExtArrayArray = ^TSingleExtArrayArray;
              TSingleExtArrayArray =  array[0..0] of PSingleExtArray;

         type PSingleExtArrayArrayArray = ^TSingleExtArrayArrayArray;
              TSingleExtArrayArrayArray =  array[0..0] of PSingleExtArrayArray;

         type PDoubleExtArray  = ^TDoubleExtArray;
              TDoubleExtArray  =  array[0..0] of double;

         type PDoubleExtArrayArray  = ^TDoubleExtArrayArray;
              TDoubleExtArrayArray  =  array[0..0] of PDoubleExtArray;

         type FileOfSingle = File of Single;
              PFileOfSingleExtArray = ^TFileOfSingleExtArray;
              TFileOfSingleExtArray =  array[0..0] of FileOfSingle;

         type PFileExtArray = ^TFileExtArray;
              TFileExtArray =  array[0..0] of File;

              //Array of String[8]
         type String8  = String[8];
              PString8ExtArray = ^TString8ExtArray;
              TString8ExtArray =  array[0..0] of String8;

              //Array of String[16]
         type String16  = String[16];
              PString16ExtArray = ^TString16ExtArray;
              TString16ExtArray =  array[0..0] of String16;

              //Array of String[32]
         type String32  = String[32];
              PString32ExtArray = ^TString32ExtArray;
              TString32ExtArray =  array[0..0] of String32;

              //Array of String[80]
         type String80  = String[80];
              PString80ExtArray = ^TString80ExtArray;
              TString80ExtArray =  array[0..0] of String80;

              //Pointer to a character
         type PAChar = ^Char;

              //Pointer to a byte
         type PByte = ^Byte;

         type T24Bit = record
                           Byte0,Byte1,Byte2:byte;
                           end;
              P24BitExtArray = ^T24BitExtArray;
              T24BitExtArray =  array[0..0] of T24Bit;

         type TComplexSingle = record
                                       R,I:single;
                                 end;

              PComplexSingleExtArray = ^TComplexSingleExtArray;
              TComplexSingleExtArray = array[0..0] of TComplexSingle;
              PComplexSingleExtArrayArray = ^TComplexSingleExtArrayArray;
              TComplexSingleExtArrayArray = array[0..0] of PComplexSingleExtArray;
              PComplexSingleExtArrayArrayArray = ^TComplexSingleExtArrayArrayArray;
              TComplexSingleExtArrayArrayArray = array[0..0] of PComplexSingleExtArrayArray;
              PComplexSingleExtArrayArrayArrayArray = ^TComplexSingleExtArrayArrayArrayArray;
              TComplexSingleExtArrayArrayArrayArray = array[0..0] of PComplexSingleExtArrayArrayArray;

         const ComplexSingleZero:TComplexSingle = (R:0;I:0);

         type TSinglePair = record
                                  a,b:single;
              end;

              TIntegerPair = record
                                  a,b:integer;
              end;

implementation
end.
