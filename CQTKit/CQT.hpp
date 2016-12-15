//
//  CQT.hpp
//  CQTKit
//
//  Created by Matt on 12/14/16.
//  Copyright © 2016 MattRajca. All rights reserved.
//

#ifndef CQT_hpp
#define CQT_hpp

namespace CQTKit {

enum class WindowFunction {
	Hamming
};

class CQT {
public:
	/// Sets up a CQT with a min/max frequency, a fixed number of bins (usually 12),
	/// a given sampling frequency, and a windowing function (currently only Hamming is supported).
	CQT(float minFreq, float maxFreq, int bins, float sampleFreq, WindowFunction);

	/// Runs the CQT forward, taking vector `x` of length `N` from the time domain to the
	/// frequency domain.
	///
	/// The returned vector is of length `k()` and must be freed by the caller.
	float *forward(float *x, int N) const;

	inline int k() const {
		return m_K;
	}

	virtual ~CQT();

private:
	float m_minFreq;
	float m_maxFreq;
	int m_bins;
	float m_sampleFreq;
	WindowFunction m_windowFunction;

	float m_Q;
	void *m_fftSetup;
	int m_fftLength;
	int m_fftLogLength;
	float *m_fftLengthMatrixReal;
	float *m_fftLengthMatrixImag;
	float *m_kernelReal;
	float *m_kernelImag;
	int m_K;
	float *m_KVectorReal;
	float *m_KVectorImag;
};

} // namespace CQTKit

#endif /* CQT_hpp */
