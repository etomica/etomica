package etomica.modules.pistoncylinder;

import etomica.Phase;
import etomica.data.meter.MeterDensity;
import etomica.potential.P1HardMovingBoundary;

/**
 * Specialized density meter for piston/cylinder system
 */
public class MeterPistonDensity extends MeterDensity {
    public MeterPistonDensity(P1HardMovingBoundary potential, int wallDim) {
        super();
        pistonPotential = potential;
        wallD = wallDim;
    }

    public double getDataAsScalar(Phase p) {
        double totDensity = super.getDataAsScalar(p);
        double d = p.boundary().dimensions().x(wallD);
        return totDensity * (d / (d - pistonPotential.getWallPosition()));
    }

    private P1HardMovingBoundary pistonPotential;
    private int wallD;
}