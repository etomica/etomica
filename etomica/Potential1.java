package etomica; 

/**
 * Potential acting on a single atom or atom group.
 *
 * @author David Kofke
 */
public abstract class Potential1 extends Potential {
      
	protected Space.Boundary boundary;
	
    public Potential1(SimulationElement parent) {
        super(1, parent);
    }

    public void setPhase(Phase phase) {
    	boundary = phase.boundary();
    }
    
    /**
     * Marker interface indicating that a one-body potential is an intramolecular
     * potential, and not, e.g., a potential of interaction with an external field.
     * This is useful when computing energy changes for molecule translations and
     * rotations, for which intramolecular contributions can be ignored.
     */
    public interface Intramolecular {}
    
	/**
	 * Methods needed to describe the behavior of a hard one-body potential.  
	 * A hard potential describes impulsive interactions, in which the energy undergoes a step
	 * change at some point in the space.
	 */
 
	 /* History of changes
	  * 01/26/03 (DAK) changed to interface and put inside Potential1 class
	  * 08/26/02 (DAK) Modified calculate to handle cases of atomCount 0 or 1 in
	  * iterator directive
	  */

}//end of Potential1



