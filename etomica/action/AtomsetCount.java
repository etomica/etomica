package etomica.action;

import etomica.api.IAtomSet;


/**
 * Action that simply counts the number of times the actionPerformed
 * method is invoked.
 */

 //used by some iterators to implement their size() method

public class AtomsetCount extends AtomsetActionAdapter {


    /**
     * Increments the call-counter by 1.
     */
    public void actionPerformed(IAtomSet atom) {callCount++;}

    /**
	 * Sets the callCount to zero.
	 */
	public void reset() {callCount = 0;}
	
	/**
	 * @return the number of times actionPerformed has been invoked
	 * since instantiation or last reset() call.
	 */
	public int callCount() {return callCount;}

    private static final long serialVersionUID = 1L;
	private int callCount = 0;
}
