//
//  CQT.cpp
//  CQTKit
//
//  Created by Matt on 12/14/16.
//  Copyright © 2016 MattRajca. All rights reserved.
//

#include "CQT.hpp"

#include <Accelerate/Accelerate.h>

namespace CQTKit {

CQT::CQT(float minFreq, float maxFreq, int bins, float sampleFreq, WindowFunction func) : m_minFreq(minFreq), m_maxFreq(maxFreq), m_bins(bins), m_sampleFreq(sampleFreq), m_windowFunction(func)
{
	const int K = static_cast<int>(ceilf(bins * log2f(maxFreq / minFreq)));
	m_K = K;
	m_Q = 1 / (powf(2, 1.0f / bins) - 1);
	m_fftLogLength = static_cast<int>(ceilf(log2f(m_Q * sampleFreq / minFreq)));
	m_fftLength = static_cast<int>(powf(2, m_fftLogLength));

	DSPSplitComplex kernel;
	kernel.realp = new float[K * m_fftLength];
	kernel.imagp = new float[K * m_fftLength];
	bzero(kernel.realp, sizeof(float) * K * m_fftLength);
	bzero(kernel.imagp, sizeof(float) * K * m_fftLength);
	m_kernelReal = kernel.realp;
	m_kernelImag = kernel.imagp;

	FFTSetup setup = vDSP_create_fftsetup(m_fftLogLength, FFT_RADIX2);
	m_fftSetup = reinterpret_cast<void *>(setup);

	const int maxN = static_cast<int>(ceilf(m_Q * m_sampleFreq / minFreq));

	DSPSplitComplex window;
	window.realp = new float[maxN];
	window.imagp = new float[maxN];

	DSPSplitComplex exps;
	exps.realp = new float[maxN];
	exps.imagp = new float[maxN];

	for (int k = K; k > 0; --k) {
		const int N = static_cast<int>(ceilf(m_Q * m_sampleFreq / (minFreq * pow(2, static_cast<float>(k - 1) / bins))));
		assert(N <= maxN);

		bzero(window.imagp, sizeof(float) * N);

		switch (m_windowFunction) {
		case WindowFunction::Hamming:
			vDSP_hamm_window(window.realp, N, 0);
			break;
		case WindowFunction::Hann:
			vDSP_hann_window(window.realp, N, 0);
			break;
		}

		for (int i = 0; i < N; ++i) {
			float inner = 2 * M_PI * m_Q * static_cast<float>(i) / N;
			exps.realp[i] = inner;
			exps.imagp[i] = inner;
		}

		vvcosf(exps.realp, exps.realp, &N);
		vvsinf(exps.imagp, exps.imagp, &N);

		const float nFloat = static_cast<float>(N);
		vDSP_vsdiv(exps.realp, 1, &nFloat, exps.realp, 1, N);
		vDSP_vsdiv(exps.imagp, 1, &nFloat, exps.imagp, 1, N);

		DSPSplitComplex complexPointer;
		complexPointer.realp = &kernel.realp[(k-1) * m_fftLength];
		complexPointer.imagp = &kernel.imagp[(k-1) * m_fftLength];

		vDSP_zvmul(&window, 1, &exps, 1, &complexPointer, 1, N < m_fftLength ? N : m_fftLength, 1);
		vDSP_fft_zip(setup, &complexPointer, 1, m_fftLogLength, FFT_FORWARD);
	}

	delete[] window.realp;
	delete[] window.imagp;

	delete[] exps.realp;
	delete[] exps.imagp;

	vDSP_mtrans(kernel.realp, 1, kernel.realp, 1, m_fftLength, K);
	vDSP_mtrans(kernel.imagp, 1, kernel.imagp, 1, m_fftLength, K);

	vDSP_zvconj(&kernel, 1, &kernel, 1, m_fftLength * K);

	const float lengthFloat = static_cast<float>(m_fftLength);
	vDSP_vsdiv(kernel.realp, 1, &lengthFloat, kernel.realp, 1, m_fftLength * K);
	vDSP_vsdiv(kernel.imagp, 1, &lengthFloat, kernel.imagp, 1, m_fftLength * K);

	m_fftLengthMatrixReal = new float[m_fftLength];
	m_fftLengthMatrixImag = new float[m_fftLength];

	m_KVectorReal = new float[m_K];
	m_KVectorImag = new float[m_K];
}

CQT::~CQT()
{
	delete[] m_fftLengthMatrixReal;
	delete[] m_fftLengthMatrixImag;

	delete[] m_KVectorReal;
	delete[] m_KVectorImag;

	vDSP_destroy_fftsetup(reinterpret_cast<FFTSetup>(m_fftSetup));
}

float *CQT::forward(float *x, int length) const {
	float *mags = new float[m_K];
	forward(x, mags, length);
	return mags;
}

void CQT::forward(const float *x, float *mags, int length) const
{
	DSPSplitComplex xMatrix;
	xMatrix.realp = m_fftLengthMatrixReal;
	xMatrix.imagp = m_fftLengthMatrixImag;

	if (length < m_fftLength) {
		bzero(xMatrix.realp, sizeof(float) * m_fftLength);
	}

	memcpy(xMatrix.realp, x, sizeof(float) * (length > m_fftLength ? m_fftLength : length));
	bzero(xMatrix.imagp, sizeof(float) * m_fftLength);

	vDSP_fft_zip(reinterpret_cast<FFTSetup>(m_fftSetup), &xMatrix, 1, m_fftLogLength, FFT_FORWARD);

	DSPSplitComplex result;
	result.realp = m_KVectorReal;
	result.imagp = m_KVectorImag;

	DSPSplitComplex existing;
	existing.realp = m_kernelReal;
	existing.imagp = m_kernelImag;
	vDSP_zmmul(&xMatrix, 1, &existing, 1, &result, 1, 1, m_K, m_fftLength);
	vDSP_zvmags(&result, 1, mags, 1, m_K);
}

} // namespace CQTKit
