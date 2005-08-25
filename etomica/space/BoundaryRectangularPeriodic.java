package etomica.space;

import etomica.simulation.Simulation;



/**
 * Rectangular boundary that is periodic in every dimension.
 */

/*
 * History
 * Created on Apr 7, 2005 by kofke
 */
public class BoundaryRectangularPeriodic extends BoundaryRectangular {

    public BoundaryRectangularPeriodic(Simulation sim) {
        this(sim.space, sim.getDefaults().boxSize);
    }
    
    /**
     * @param space
     */
    public BoundaryRectangularPeriodic(Space space, double boxSize) {
        super(space, makePeriodicity(space.D()), boxSize);
    }

    public void nearestImage(Vector dr) {
        dr.PE(dimensionsHalf);
        dr.mod(dimensions);
        dr.ME(dimensionsHalf);
    }

    public Vector centralImage(Vector r) {
        modShift.EModShift(r, dimensions);
        return modShift;
    }

    private static boolean[] makePeriodicity(int D) {
        boolean[] isPeriodic = new boolean[D];
        for (int i=0; i<D; i++) {
            isPeriodic[i] = true;
        }
        return isPeriodic;
    }
    
}
