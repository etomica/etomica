package etomica.space;

import etomica.simulation.Simulation;


/**
 * Class for implementing slit periodic boundary conditions, in which
 * one dimension is not periodic.
 */

/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class BoundaryRectangularSlit extends BoundaryRectangular {
    
    public BoundaryRectangularSlit(Simulation sim) {
        //consumer can set appropriate slit dim later
        this(sim,0);
    }
    
    public BoundaryRectangularSlit(Simulation sim, int slitDim) {
        this(sim.space, slitDim, sim.getDefaults().boxSize);
    }
    
    /**
     * Constructor for periodic boundary conditions with a slit 
     * in the given dimension.
     * @param space
     * @param slitDim slit dimension (in which PBC are not imposed).
     */
    public BoundaryRectangularSlit(Space space, int slitDim, double boxSize) {
        super(space,makePeriodicity(space.D(),slitDim),boxSize);
        sDim = slitDim;
    }
    
    public void setSlitDim(int slitDim) {
        isPeriodic = makePeriodicity(space.D(),slitDim);
        sDim = slitDim;
    }
    
    public int getSlitDim() {
        return sDim;
    }
    
    public void nearestImage(Vector dr) {
        double x = dr.x(sDim);
        dr.PE(dimensionsHalf);
        dr.mod(dimensions);
        dr.ME(dimensionsHalf);
        dr.setX(sDim,x);
    }
    public Vector centralImage(Vector r) {
        //superclass sets modShift
        modShift.EModShift(r, dimensions);
        modShift.setX(sDim,0.0);
        return modShift;
    }

    private static boolean[] makePeriodicity(int D, int slitDim) {
        boolean[] isPeriodic = new boolean[D];
        for (int i=0; i<D; i++) {
            isPeriodic[i] = slitDim != i;
        }
        return isPeriodic;
    }
    
    private int sDim;
}
