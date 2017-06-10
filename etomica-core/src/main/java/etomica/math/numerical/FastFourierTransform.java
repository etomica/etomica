/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.math.numerical;

/**
 *
 * This utility receives a set of data points, either Real, Imaginary or Both
 * and performs a Discrete Fourier Transform on them using the Fast Fourier
 * Transform algorithm. The class modifies the passed arrays directly, however
 * they can use replicate copies instead using the BACKUP toggle. 
 * 
 * The class also maintains an array of the frequency indexes of the respective 
 * datum. This index can be used for graphing purposes (it would be the x axis, 
 * while the y axis is the real and/or imaginary transformed values.
 *  
 * @author Mitchell Tai
 */

/*
 * Created on Jun 25, 2003
 */
public class FastFourierTransform implements java.io.Serializable {

	private double[] real;					// array of real plots
	private double[] imaginary;				// array of imaginary plots
	public int length;						// length of array
	
	/**
	 * Constructors accepting: nothing; all reals; both real and imaginary
	 */
	public FastFourierTransform() {
	}
	public FastFourierTransform(double[] dataReal) {
		setData(dataReal);
	}	
	public FastFourierTransform(double[] dataReal, double[] dataImaginary) {
		setData(dataReal,dataImaginary);
	}
	
	/**
	 * Saves the two passed arrays of doubles into this class, real and 
	 * imaginary in order to transform. Performs additional checks and indexes 
	 * the array.
	 */
	public void setData(double[] dataReal, double[] dataImaginary) {
	 // Check to see if both arrays match in length
		if(dataReal.length != dataImaginary.length)
			throw new IllegalArgumentException ("Passed arrays not of equal size");
		
		length = dataReal.length;
		
	 // Check to see if value n is a power of 2 (log2 of n is a whole number)
		if ((Math.log(length)/Math.log(2))%1 != 0) {
			throw new IllegalArgumentException("Array Index is not a power of 2");
		}
		
		real=dataReal;
		imaginary=dataImaginary;
	}

	/**
	 * For passing only all reals (true)  or all imaginaries (false)
	 */
	public void setData(double[] data, boolean real){
		double[] emptySet = new double [data.length];
		for (int i=0;i<emptySet.length;i++) emptySet[i]=0;

		if(real) setData(data, emptySet);
		else setData(emptySet,data);
	}	

	/**
	 * For passing only one set of data - assumes passing all reals
	 */
	public void setData(double[] data){setData(data, true);}

	/**
	 * General Transform function that runs both the forward and reverse 
	 * Transform by passing +1 or -1 isign.
	 * Works at Nlog2N speed
	 */
	private void FFT(int isign) {

		int j,k,i2;
  		double temp;

		int nn=real.length;
		int m = (int)(Math.log(nn)/Math.log(2));

   		/* Do the bit reversal */
 		i2 = nn >> 1;
		j = 0;
		for (int i=0;i<nn-1;i++) {
			if (i < j) {
				temp = real[i];
				real[i] = real[j];
				real[j] = temp;

				temp = imaginary[i];
				imaginary[i] = imaginary[j];
				imaginary[j] = temp;
			}
			k = i2;
			while (k <= j) {
				j -= k;
				k >>= 1;
			}
				j += k;
		}
		
		/* Compute the FFT */
		double wr,wi,wpr,wpi,wtemp,tempr,tempi;
		int mmax,istep,i1;

		wpr = -1.0; 
		wpi = 0.0;
		istep = 1;
		for (int l=0;l<m;l++) {
			mmax = istep;
			istep <<= 1;
			wr = 1.0; 
			wi = 0.0;
			for (j=0;j<mmax;j++) {
				for (int i=j;i<nn;i+=istep) {
					i1 = i + mmax;

					tempr = wr * real[i1] - wi * imaginary[i1];
					tempi = wr * imaginary[i1] + wi * real[i1];
					real[i1] = real[i] - tempr; 
					imaginary[i1] = imaginary[i] - tempi;
					real[i] += tempr;
					imaginary[i] += tempi;
				}
				wtemp =  wr * wpr - wi * wpi;
				wi = wr * wpi + wi * wpr;
				wr = wtemp;
			}
			wpi = Math.sqrt((1.0 - wpr) / 2.0);
			if (isign == 1) wpi = -wpi;
			wpr = Math.sqrt((1.0 + wpr) / 2.0);
		}
		/* Scaling for forward transform */
		if (isign == 1) {
			for (int i=0;i<nn;i++) {
				real[i] /= nn;
				imaginary[i] /= nn;
			}
		}
	}

	/**
	 * Forward Transformation of discrete data sets
	 */
	public void transform() { FFT(1);}

	/**
	 * Reverse Transformation of discrete data sets
	 */
	public void invert() { FFT(-1);}
	
	/**
	 * array value calls
	 */
	public double[] getReal() {return real;}
	public double[] getImaginary() {return imaginary;}

// in order to test the fft....	
	public static void main (String [] arg) {
	
		int samplingRate=16; 
		double[] data1 = new double [samplingRate];
		double[] data2 = new double [samplingRate];
		for(int i=0;i<data1.length;i++) { 
			data1[i] = Math.sin(Math.PI*(i+7)/8);
			data2[i] = data1[i];
		}
		System.out.println("\n");
	    for (int i=0;i<data1.length;i++) {
	   		System.out.print("[" + round(data1[i]) + ",");
	   		System.out.print(round(data2[i]) + "]\n");
		}
		
		FastFourierTransform fourier = new FastFourierTransform();
		
		fourier.setData(data1,data2);
		fourier.transform();
		//fourier.invert();

		System.out.println("\n");
	    for (int i=0;i<fourier.length;i++) {
	   		System.out.print("(" + round(fourier.getReal()[i]) + ",");
	   		System.out.print(round(fourier.getImaginary()[i]) + ":");
	   		System.out.print("\n");
	    }
	}

	// mini rounding function - only necessary for display purposes
	public static double round(double x) {
		x=Math.round(x*1000);
		return x/1000;
	}
	
}
