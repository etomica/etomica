package etomica.utility;

/**
 * Function interface that takes an array of doubles and returns a 
 * transformed array of doubles
 */

public interface Transform {

	public double[] f(double[] x);
	public double[] xHat();

	public static class Fourier implements Transform {
		FastFourierTransform fourier = new FastFourierTransform();
		public static boolean REAL=true;
		public double[] f(double[] x) {
			fourier.setData(x);
			fourier.transform();
			if (REAL) return fourier.getReal();
			return fourier.getImaginary();
		}
		public double[] xHat() {
			return fourier.getIndex();
		}
	}
}
