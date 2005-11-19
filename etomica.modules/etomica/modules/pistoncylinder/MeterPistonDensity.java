package etomica.modules.pistoncylinder;

import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.phase.Phase;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.DimensionRatio;

/**
 * Specialized density meter for piston/cylinder system
 */
public class MeterPistonDensity extends DataSourceScalar implements Meter {
    public MeterPistonDensity(P1HardMovingBoundary potential, int wallDim, double atomDiameter) {
        super("Density",new DimensionRatio(Dimension.QUANTITY,Dimension.VOLUME));
        pistonPotential = potential;
        wallD = wallDim;
        collisionDiameter = atomDiameter;
    }

    public void setAtomDiameter(double d) {
        collisionDiameter = d;
    }
    
    public double getDataAsScalar() {
        double volume = 1;
        final Vector dimensions = phase.getBoundary().getDimensions();
        for (int i=0; i<dimensions.D(); i++) {
            double d = dimensions.x(i) - collisionDiameter;
            if (i==wallD) {
                d -= pistonPotential.getWallPosition();
            }
            volume *= d;
        }
        return phase.moleculeCount()/volume;
    }
    
    public Dimension getDimension() {
        return new DimensionRatio(Dimension.QUANTITY, Dimension.VOLUME);
    }

    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
    private P1HardMovingBoundary pistonPotential;
    private int wallD;
    private double collisionDiameter;
}