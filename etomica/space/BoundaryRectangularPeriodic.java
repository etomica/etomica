package etomica.space;

import etomica.simulation.Simulation;

/**
 * Rectangular boundary that is periodic in every dimension.
 */
public class BoundaryRectangularPeriodic extends BoundaryRectangular {

    /**
     * Constructs cubic boundary with the default box-size given by the Simulation.
     */
    public BoundaryRectangularPeriodic(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().boxSize);
    }
    
    /**
     * Constructs cubic boundary for the given Space, with each edge of length boxSize.
     */
    public BoundaryRectangularPeriodic(Space space, double boxSize) {
        super(space, makePeriodicity(space.D()), boxSize);
        dimensionsHalf = space.makeVector();
        tempImage = (IVectorRandom)space.makeVector();
        // call updateDimensions again so dimensionsHalf is updated
        updateDimensions();
    }

    public void updateDimensions() {
        super.updateDimensions();
        // superclass constructor calls this before dimensionsHalf has been instantiated
        if (dimensionsHalf != null) {
            dimensionsHalf.Ea1Tv1(0.5,dimensions);
        }
    }
    
    public void nearestImage(IVector dr) {
        dr.PE(dimensionsHalf);
        dr.mod(dimensions);
        dr.ME(dimensionsHalf);
    }

    public IVector centralImage(IVector r) {
        tempImage.E(r);
        nearestImage(tempImage);
        tempImage.ME(r);
        return tempImage;
    }

    private static boolean[] makePeriodicity(int D) {
        boolean[] isPeriodic = new boolean[D];
        for (int i=0; i<D; i++) {
            isPeriodic[i] = true;
        }
        return isPeriodic;
    }
    
    private static final long serialVersionUID = 1L;
    protected final IVector dimensionsHalf;
    protected final IVectorRandom tempImage;
}
