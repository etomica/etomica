package etomica.modules.pistoncylinder;

import etomica.atom.AtomSet;
import etomica.atom.AtomsetFilter;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.ApiFiltered;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.phase.Phase;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.Boundary;
import etomica.space.IVector;

/**
 * Our own ApiFiltered that's phase-dependent
 */
public class ApiFilteredCylinder extends ApiFiltered implements AtomsetIteratorPhaseDependent {
    public ApiFilteredCylinder(AtomsetIterator iterator, AtomsetFilter filter) {
        super(iterator, filter);
    }

    public void setPhase(Phase newPhase) {
        ((AtomsetIteratorPhaseDependent)iterator).setPhase(newPhase);
    }
    private static final long serialVersionUID = 1L;

    /**
     * Filter to expclude any pair with an atom within some distance from a 
     * wall. 
     */
    public static class AtomFilterInCylinder implements AtomsetFilter {
        public AtomFilterInCylinder(Boundary boundary, P1HardMovingBoundary pistonPotential, double padding) {
            dimensions = boundary.getDimensions();
            this.pistonPotential = pistonPotential;
            this.padding = padding;
            // bit flipper goes back and forth between 1 and 2
            bitFlipper = 1;
        }
        
        public boolean accept(AtomSet atoms) {
            double radius = pistonPotential.getCollisionRadius()+padding;
            // always reject if both atoms are near a wall.  always accept if
            // both atoms are away from the wall.  If one is near and one not, 
            // accept the pair half the time.  RDF needs this to avoid 
            // over-counting pairs with one near the wall.  Ideally, we'd 
            // accept them all and weight them half as much. 
            int numOut = 0;
            for (int i=0; i<2; i++) {
                IVector pos = ((IAtomPositioned)atoms.getAtom(i)).getPosition();
                if (pos.x(0) < -0.5*dimensions.x(0)+radius ||
                    pos.x(0) >  0.5*dimensions.x(0)-radius ||
                    pos.x(1) >  0.5*dimensions.x(1)-radius ||
                    pos.x(1) < pistonPotential.getWallPosition()+radius) {
                    numOut++;
                }
            }
            // twiddle the last two bits, 1=>2, 2=>1
            // numOut=0 is always accepted, numOut=2 is never accepted
            bitFlipper ^= 3;
            return numOut < bitFlipper;
        }
        
        private double padding;
        private final IVector dimensions;
        private final P1HardMovingBoundary pistonPotential;
        private int bitFlipper;
    }
}