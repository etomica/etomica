/*
 * History
 * Created on Jul 28, 2004 by kofke
 */
package etomica;

/**
 * @author kofke
 *
 * Defines a conversion between an array of double and an Object.
 * The double[] form is used internally for accumulating averages, and
 * the translator converts this simple data structure into an object
 * (such as a Space.Vector) that would more appropriately represent
 * the data.
 */
public interface DataTranslator {

	/**
	 * Converts the given array of doubles into the object.
	 */
	public Object fromArray(double[] x);
	
	/**
	 * Converts the object into an array of doubles.
	 */
	//perhaps we don't want this method as part of interface
	public double[] toArray(Object obj);

	
	/**
	 * Simple identity translator, in which the conversion Object is the 
	 * double[] array itself.
	 */
	public static final DataTranslator IDENTITY = new DataTranslator() {
		public final Object fromArray(double[] x) {return x;}
		public final double[] toArray(Object obj) {return (double[])obj;}
	};
}
