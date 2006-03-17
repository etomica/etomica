package etomica.modules.pistoncylinder;

import etomica.data.DataSourceScalar;
import etomica.data.meter.Meter;
import etomica.phase.Phase;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.Vector;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;

/**
 * Specialized density meter for piston/cylinder system
 */
public class MeterPistonDensity extends DataSourceScalar implements Meter {
    public MeterPistonDensity(P1HardMovingBoundary potential, int wallDim, double atomDiameter) {
        super("Density",new DimensionRatio(Quantity.DIMENSION, Volume.DIMENSION));
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
            double d = dimensions.x(i)*0.5 - collisionDiameter;
            if (i==wallD) {
                d -= pistonPotential.getWallPosition();
            }
            else {
            		d += dimensions.x(i)*0.5;
            }
            volume *= d;
        }
        return phase.moleculeCount()/volume;
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