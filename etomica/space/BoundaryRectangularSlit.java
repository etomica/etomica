package etomica.space;

import etomica.simulation.Simulation;
import etomica.util.IRandom;

/**
 * Class for implementing slit periodic boundary conditions, in which
 * one dimension is not periodic.  Selection of non-periodic dimension may be changed
 * after construction.
 */
public class BoundaryRectangularSlit extends BoundaryRectangular {
    
    /**
     * Makes cubic volume with the x-dimension (index 0)
     * not periodic.  Length of each box edge is given by default boxSize in
     * given Simulation.
     */
    public BoundaryRectangularSlit(Simulation sim) {
        //consumer can set appropriate slit dim later
        this(sim,0);
    }
    
    /**
     * Makes cubic volume with the indicated dimension not periodic.
     * Length of each box edge is given by the default boxSize in the
     * given Simulation.
     * 
     * @param slitDim index indicating dimension that is not periodic (0 for x-dimension,
     * 1 for y-dimension, etc.).
     * @throws IllegalArgumentException if not (0 <= slitDim < space.D).
     */
    public BoundaryRectangularSlit(Simulation sim, int slitDim) {
        this(sim.getSpace(), sim.getRandom(), slitDim, sim.getDefaults().boxSize);
    }
    
    /**
     * Constructor for periodic boundary conditions with a slit 
     * in the given dimension.
     * @param space
     * @param slitDim slit dimension (in which PBC is not imposed).
     */
    public BoundaryRectangularSlit(Space space, IRandom random, int slitDim, double boxSize) {
        super(space,random,makePeriodicity(space.D(),slitDim),boxSize);
        sDim = slitDim;
        dimensionsHalf = space.makeVector();
        tempImage = space.makeVector();
        // call updateDimensions again so dimensionsHalf is updated
        updateDimensions();
    }

    /**
     * Sets the non-periodic dimension to that indicated by the given value 
     * (0 is x-dimension, 1 is y-dimension, etc.).
     * 
     * @throws IllegalArgumentException if not (0 <= slitDim < space.D).
     */
    public void setSlitDim(int slitDim) {
        isPeriodic = makePeriodicity(space.D(),slitDim);
        sDim = slitDim;
    }
    
    /**
     * Returns index of non-periodic dimension.
     */
    public int getSlitDim() {
        return sDim;
    }
    
    protected void updateDimensions() {
        super.updateDimensions();
        // superclass constructor calls this before dimensionsHalf has been instantiated
        if (dimensionsHalf != null) {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
        }
    }

    public void nearestImage(IVector dr) {
        double x = dr.x(sDim);
        dr.PE(dimensionsHalf);
        dr.mod(dimensions);
        dr.ME(dimensionsHalf);
        dr.setX(sDim,x);
    }
    
    public IVector centralImage(IVector r) {
        tempImage.E(r);
        nearestImage(tempImage);
        tempImage.ME(r);
        return tempImage;
    }

    private static boolean[] makePeriodicity(int D, int slitDim) {
        if(slitDim >= D || slitDim < 0) { 
            throw new IllegalArgumentException("Indicated dimension ("+slitDim+") is not valid for the space");
        }
        boolean[] isPeriodic = new boolean[D];
        for (int i=0; i<D; i++) {
            isPeriodic[i] = slitDim != i;
        }
        return isPeriodic;
    }
    
    private int sDim;
    private static final long serialVersionUID = 1L;
    protected final IVector dimensionsHalf;
    protected final IVector tempImage;
}
