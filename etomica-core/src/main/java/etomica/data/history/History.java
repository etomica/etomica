/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * History
 * Created on Aug 9, 2004 by kofke
 */
package etomica.data.history;


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
     */
	public int getHistoryLength();
	
	/**
	 * Returns the number of data samples currently held by the History.
	 */
	public int getSampleCount();

    /**
     * Resets the History.  Memory of all data given via addValue is dropped.
     */
	public void reset();

    /**
     * Adds the given x, y pair to the History.  Returns true if getHistory
     * will now return different values.
     */
    public boolean addValue(double x, double y);

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
}
