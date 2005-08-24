package etomica.potential; 

import etomica.Phase;
import etomica.space.Boundary;
import etomica.space.Space;

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
      
	protected Boundary boundary;
	
    public Potential1(Space space) {
        super(1, space);
    }

    public void setPhase(Phase phase) {
    	boundary = phase.boundary();
    }
    
    /**
     * Returns zero.
     */
    public double getRange() {
        return 0.0;
    }
    
    /**
     * Marker interface indicating that a one-body potential is an intramolecular
     * potential, and not, e.g., a potential of interaction with an external field.
     * This is useful when computing energy changes for molecule translations and
     * rotations, for which intramolecular contributions can be ignored.
     */
    public interface Intramolecular {}
     
}//end of Potential1



