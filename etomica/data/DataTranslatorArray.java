/*
 * Created on Jul 29, 2004
 *
 * 
 * 
 */
package etomica.data;

import etomica.DataTranslator;

/**
 * @author nancycribbin
 *
 */
public class DataTranslatorArray implements DataTranslator {
	
	/**
	 * For the 2 dimensional array, the first dimension (number of rows) is dimension1.  
	 * The second dimension (number of columns) is dimension2.
	 */
	private int dimension1, dimension2;
	
	/**
	 * Dimension is the length of the 1 dimensional array.
	 */
	private int dimension;
	
	private double[][] array2D;
	private double[] array1D;
	
	/**
	 * Constructor.
	 */
	public DataTranslatorArray(int dimension1, int dimension2) {
		setDimensions(dimension1, dimension2);
	}
	/**
	 * Converts a one dimensional array into a two dimensional array.
	 * @param x
	 * @return
	 */
	public Object fromArray(double[] x){
		if(array2D == null) {array2D = new double[dimension1][dimension2];}
		
		int k = 0;
		for(int i = 0; i < dimension1; i++) {
			System.arraycopy(x, k, array2D[i], 0, dimension2);
			k += dimension2;
//			for(int j = 0; j < dimension2; j++){
//				array2D[i][j] = x[k++];
//			}
		}
		return array2D;
	}
	
	/**
	 * Converts a two dimensional array into a one dimensional array.
	 * If 2D array has only one row, returns that row as the 1D array; 
	 * otherwise copies each row consecutively into a new (at first call)
	 * 1D array.
	 * @param obj
	 * @return
	 */
	public double[] toArray(Object obj){
		
		if(dimension1 == 1) return ((double[][])obj)[0];
		
		if(array1D == null) {array1D = new double[dimension];}
		
		double[][] castObject = (double[][]) obj;
		
		int k = 0;
		for(int i = 0; i < dimension1; i++) {
			System.arraycopy(castObject[i], 0, array1D, k, dimension2);
			k += dimension2;
//			for(int j = 0; j < dimension2; j++){
//				array1D[k++] = castObject[i][j];
//			}
		}
		return array1D;
	}
	

	public static void main(String[] args) {
//		DataTranslatorArray test = new DataTranslatorArray();
		double[][] door = new double[3][4];
//		test.toArray(door);
		
	}

	/**
	 * @param dimension1 The dimension1 to set.
	 * @param dimension2 The dimension2 to set.
	 */
	public void setDimensions(int dimension1, int dimension2) {
		this.dimension1 = dimension1;
		this.dimension2 = dimension2;
		this.dimension = dimension1 * dimension2;
		
		array1D = null;
		array2D = null;
	}
	/**
	 * @return Returns the dimension1.
	 */
	public int getDimension1() {
		return dimension1;
	}
	/**
	 * @return Returns the dimension2.
	 */
	public int getDimension2() {
		return dimension2;
	}

}
