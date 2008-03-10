package etomica.modules.pistoncylinder;

import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardMoleculeMonatomic;
import etomica.potential.P1HardMovingBoundary;
import etomica.space.Space;

/**
 * data source front for virial sum from P1HardMovingBoundary
 * returns pressure exerted on the wall by atoms
 */
public class DataSourceWallPressure extends MeterPressureHard {

    public DataSourceWallPressure(Space space, P1HardMoleculeMonatomic potentialWrapper) {
        super(space);
        wallPotentialWrapper = potentialWrapper;
    }
    
    /**
     * Implementation of CollisionListener interface
     * Adds collision virial (from potential) to accumulator
     */
    public void collisionAction(IntegratorHard.Agent agent) {
        if (agent.collisionPotential == wallPotentialWrapper) {
            if (wallPotentialWrapper.getSpace().D() == 2) {
                virialSum += ((P1HardMovingBoundary)wallPotentialWrapper.getWrappedPotential()).lastWallVirial();
            }
            else {
                virialSum -= ((P1HardMovingBoundary)wallPotentialWrapper.getWrappedPotential()).lastWallVirial();
            }
        }
    }
    
    public double getDataAsScalar() {
        double currentTime = integratorHard.getCurrentTime();
        double value = virialSum / (currentTime - lastTime);
        lastTime = currentTime;
        virialSum = 0;
        return value;
    }
    
    private static final long serialVersionUID = 1L;
    protected final P1HardMoleculeMonatomic wallPotentialWrapper;
}
