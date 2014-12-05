/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.util;

/**
 * Function interface that takes an array of doubles and returns a 
 * transformed array of doubles
 */

public interface Transform {

	public double[] f(double[] x);

	public static class Fourier implements Transform {
		FastFourierTransform fourier = new FastFourierTransform();
		public static boolean REAL=true;
		public double[] f(double[] x) {
			fourier.setData(x);
			fourier.transform();
			if (REAL) return fourier.getReal();
			return fourier.getImaginary();
		}
	}
}
