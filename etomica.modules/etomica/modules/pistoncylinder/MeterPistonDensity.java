package etomica.modules.pistoncylinder;

import etomica.Phase;
import etomica.data.meter.MeterScalar;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;

/**
 * Specialized density meter for piston/cylinder system
 */
public class MeterPistonDensity extends MeterScalar {
    public MeterPistonDensity(P1HardMovingBoundary potential, int wallDim, double atomDiameter) {
        super();
        pistonPotential = potential;
        wallD = wallDim;
        collisionDiameter = atomDiameter;
    }

    public double getDataAsScalar(Phase p) {
        double volume = 1;
        final Vector dimensions = p.boundary().dimensions();
        int D = dimensions.D();
        for (int i=0; i<dimensions.D(); i++) {
            double d = dimensions.x(i) - collisionDiameter;
            if (i==wallD) {
                d -= pistonPotential.getWallPosition();
            }
            volume *= d;
        }
        return p.moleculeCount()/volume;
    }
    
    public Dimension getDimension() {
        return new DimensionRatio(Dimension.QUANTITY, Dimension.VOLUME);
    }

    private P1HardMovingBoundary pistonPotential;
    private int wallD;
    private double collisionDiameter;
}