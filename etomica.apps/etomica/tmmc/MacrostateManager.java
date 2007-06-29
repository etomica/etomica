package etomica.tmmc;

import etomica.box.Box;

/**
 * Interface for class that defines the macrostates.  Provides
 * methods that tell the number of expected macrostates in the
 * box, and that give an index indicating which macrostate
 * the box presently is in.
 */
public interface MacrostateManager {
 /**
  * Number of expected macrostates in the box.
  */
    public int numberOfStates(Box p);
 /**
  * Returns an index identifying the current macrostate occupied by the box.
  */
    public int stateIndex(Box p);
    
 /**
  * Returns the value of the macrostate variable corresponding to the given index.
  */
    public double state(int i);
}//end of MacrostateManager
