using System;
using System.Runtime.InteropServices;
using TINS;

namespace Test
{
	class Program
	{
		const string dllpath = @"C:\Users\hbarz\source\repos\superlets_wrapper\x64\Release\superlets_wrapper.dll";

		[DllImport(dllpath)]
		static extern unsafe IntPtr asrwt_alloc(int input_size,				float sampling_rate,		
												float freq_low,				float freq_high,
												int freq_count,				float cycle_count,		
												int superresolution_low,	int superresolution_high,	bool multiplicative);
		[DllImport(dllpath)]
		static extern unsafe int asrwt_execute(IntPtr asrwt, float* input, float* output);
		[DllImport(dllpath)]
		static extern unsafe void asrwt_free(IntPtr asrwt);


		static unsafe void Main(string[] args)
		{
			int N					= 1000;
			float FS				= 1000f;
			
			var data				= new RowVector<float>(N);
			var spectrum			= new Matrix<float>(41, N);
			data.Sub(300, 400).Set	= Signals.Sine(40, 400, FS);

			IntPtr asrwt = asrwt_alloc(N, FS, 20, 60, 41, 2f, 15, 15, true);
			fixed (float* pIn = data.GetSpan())
			fixed (float* pOut = spectrum.GetSpan())
			{
				asrwt_execute(asrwt, pIn, pOut);
			}
			asrwt_free(asrwt);

			spectrum.Plot2D();
		}
	}
}
