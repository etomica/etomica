/*
 * History
 * Created on Aug 23, 2004 by kofke
 */
package etomica.action;

import etomica.atom.AtomSet;


/**
 * Action that simply counts the number of times the actionPerformed
 * method is invoked.
 */

 //used by some iterators to implement their size() method

public class AtomsetCount extends AtomsetActionAdapter {

    /**
     * Increments the call-counter by 1.
     */
    public void actionPerformed(AtomSet atom) {callCount++;}

    /**
	 * Sets the callCount to zero.
	 */
	public void reset() {callCount = 0;}
	
	/**
	 * @return the number of times actionPerformed has been invoked
	 * since instantiation or last reset() call.
	 */
	public int callCount() {return callCount;}

	private int callCount = 0;
}
