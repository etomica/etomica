/*
 * History
 * Created on Aug 9, 2004 by kofke
 */
package etomica.utility;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public interface History {
	public void setHistoryLength(int n);

	public int getHistoryLength();

	public void reset();

	public void addValue(double x);

	public double[] xValues();

	public double[] getHistory();
	
	/**
	 * Interface for a class that can make a History instance. An instance
	 * of such a factory is needed to instantiate an AccumulatorHistory
	 * object.
	 */
	public interface Factory {
		public History makeHistory();
		public History makeHistory(int n);
	}

}