package etomica.modules.pistoncylinder;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomPairFilter;
import etomica.atom.iterator.ApiFiltered;
import etomica.atom.iterator.AtomPairIterator;
import etomica.atom.iterator.AtomsetIteratorPhaseDependent;
import etomica.data.Data;
import etomica.data.meter.MeterRDF;
import etomica.phase.Phase;
import etomica.potential.P1HardMovingBoundary;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Vector;

/**
 * Our own MeterRDF to account for the excluded volume near the wall
 * @author andrew
 */
public class MeterRDFCylinder extends MeterRDF {
    public MeterRDFCylinder(Boundary boundary, P1HardMovingBoundary pistonPotential) {
        super(boundary.getSpace());
        this.boundary = boundary;
        this.pistonPotential = pistonPotential;
    }
    
    public Data getData() {
        super.getData();
        double boundaryVolume = boundary.volume();
        Vector dimensions = boundary.getDimensions();
        double radius = pistonPotential.getCollisionRadius();
        double rdfVolume = (dimensions.x(0)-2*radius)*(0.5*dimensions.x(1)-pistonPotential.getWallPosition()-2*radius);
        data.TE(rdfVolume/boundaryVolume);
        return data;
    }

    private static final long serialVersionUID = 1L;
    private final Boundary boundary;
    private final P1HardMovingBoundary pistonPotential;

    /**
     * Filter to expclude any pair with an atom within some distance from a 
     * wall. 
     */
    public static class AtomFilterInCylinder implements AtomPairFilter {
        public AtomFilterInCylinder(Boundary boundary, P1HardMovingBoundary pistonPotential, double padding) {
            dimensions = boundary.getDimensions();
            this.pistonPotential = pistonPotential;
            this.padding = padding;
        }
        
        public boolean accept(AtomPair atoms) {
            double radius = pistonPotential.getCollisionRadius()+padding;
            // always reject if both atoms are near a wall.  always accept if
            // both atoms are away from the wall.  If one is near and one not, 
            // accept the pair half the time.  RDF needs this to avoid 
            // over-counting pairs with one near the wall.  Ideally, we'd 
            // accept them all and weight them half as much. 
            int numOut = 0;
            for (int i=0; i<2; i++) {
                Vector pos = ((AtomLeaf)atoms.getAtom(i)).coord.position();
                if (pos.x(0) < -0.5*dimensions.x(0)+radius ||
                    pos.x(0) >  0.5*dimensions.x(0)-radius) {
                    numOut++;
                    continue;
                }
                if (pos.x(1) >  0.5*dimensions.x(1)-radius) {
                    numOut++;
                    continue;
                }
                if (pos.x(1) < pistonPotential.getWallPosition()+radius) {
                    numOut++;
                    continue;
                }
            }
            return numOut < Simulation.random.nextInt(2)+1;
        }
        
        private double padding;
        private final Vector dimensions;
        private final P1HardMovingBoundary pistonPotential;
    }

    /**
     * Our own ApiFiltered that's phase-dependent
     */
    public static class ApiFilteredCylinder extends ApiFiltered implements AtomsetIteratorPhaseDependent {
        public ApiFilteredCylinder(AtomPairIterator iterator, AtomPairFilter filter) {
            super(iterator, filter);
        }

        public void setPhase(Phase newPhase) {
            ((AtomsetIteratorPhaseDependent)iterator).setPhase(newPhase);
        }
        private static final long serialVersionUID = 1L;
    }
}