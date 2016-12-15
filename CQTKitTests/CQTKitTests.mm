//
//  CQTKitTests.m
//  CQTKitTests
//
//  Created by Matt on 12/14/16.
//  Copyright © 2016 MattRajca. All rights reserved.
//

#import <XCTest/XCTest.h>

#import <CQTKit/CQT.hpp>

@interface CQTKitTests : XCTestCase
@end

@implementation CQTKitTests

using namespace CQTKit;

static constexpr int sampleRate = 44100;
static constexpr int duration = 1;

- (void)testBasics
{
	// Sets things up with a minimum frequency of 40 Hz, max of 22050 Hz, 12 bins, a sample rate
	// of 44100 Hz, and a Hamming window function.
	CQT qt(40, 22050, 12, sampleRate, WindowFunction::Hamming);

	// Generates a tone at 440 Hz.
	float *x = new float[sampleRate * duration];
	for (int i = 0; i < sampleRate * duration; i++) {
		x[i] = sinf(2 * M_PI * i * (440.0f / sampleRate));
	}

	// Runs the CQT forward (from the time to frequency domain).
	float *results = qt.forward(x, 44100);
	delete[] x;

	// Finds the most likely tone.
	float max = results[0];
	int idx = 0;
	for (int i = 1; i < qt.k(); i++) {
		if (results[i] > max) {
			max = results[i];
			idx = i;
		}
	}

	delete[] results;

	// Ensures we have the highest response at A4. The offset is due to us starting at 40 Hz.
	XCTAssertEqual(idx + 27, 69 /* A4 */);
}

@end
