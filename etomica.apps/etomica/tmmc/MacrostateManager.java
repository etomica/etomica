package etomica.tmmc;

import etomica.*;

/**
 * Interface for class that defines the macrostates.  Provides
 * methods that tell the number of expected macrostates in the
 * phase, and that give an index indicating which macrostate
 * the phase presently is in.
 */
public interface MacrostateManager {
 /**
  * Number of expected macrostates in the phase.
  */
    public int numberOfStates(Phase p);
 /**
  * Returns an index identifying the current macrostate occupied by the phase.
  */
    public int stateIndex(Phase p);
    
 /**
  * Returns the value of the macrostate variable corresponding to the given index.
  */
    public double state(int i);
}//end of MacrostateManager
