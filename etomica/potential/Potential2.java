package etomica.potential; 

import etomica.Phase;
import etomica.Potential;
import etomica.Space;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.NeighborCriterionAll;
import etomica.space.CoordinatePair;

/**
 * Potential acting on or within an atom, or between a pair of atoms or atom
 * groups. Contains other Potential instances that describe the specific
 * interactions between the atoms of the group(s).
 *
 * @author David Kofke
 */

 /* History of changes
  * 07/13/02 (DAK) Restructured instantiation of LRC potential
  * 07/15/02 (DAK) Constructor makes P0LRC only if instance of Potential2SoftSpherical
  * 12/06/02 (DAK) Added setIterators1A method
  * 01/27/03 (DAK) Numerous changes with redesign of Potential.
  * 07/17/03 (DAK) Made calculate method not final (so etomica.virial.P2Cluster
  * could override it)
  */

public abstract class Potential2 extends Potential {
  
    public Potential2(Space space) {
        super(2, space);
        cPair = space.makeCoordinatePair();
    }

    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.boundary());
    }
    
    public void setCriterion(NeighborCriterion criterion) {
        this.criterion = criterion;
    }
    public NeighborCriterion getCriterion() {
        return criterion;
    }
    
    protected final CoordinatePair cPair;
    protected NeighborCriterion criterion = new NeighborCriterionAll();
}//end of Potential2



