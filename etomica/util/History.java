/*
 * History
 * Created on Aug 9, 2004 by kofke
 */
package etomica.util;


/**
 * Interface for an object that keeps track of x and y values and return
 * their history when asked.  Behavior when the history has more data than
 * the set length length is determined by the subclasses.
 */
public interface History {

    /**
     * Sets the number of data points to be tracked by the History.
     */
	public void setHistoryLength(int n);

    /**
     * Returns the number of data points tracked by the History.
     * @return
     */
	public int getHistoryLength();

    /**
     * Resets the History.  Memory of all data given via addValue is dropped.
     */
	public void reset();

    /**
     * Adds the given x, y pair to the History.
     */
    public void addValue(double x, double y);

    /**
     * Returns an array containing the x values that have been added to the
     * History.
     */
	public double[] getXValues();

    /**
     * Returns an array containing the y values that have been added to the 
     * History.
     */
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