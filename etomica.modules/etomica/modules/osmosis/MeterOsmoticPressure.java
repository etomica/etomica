package etomica.modules.osmosis;
import etomica.EtomicaInfo;
import etomica.data.meter.MeterPressureHard;
import etomica.integrator.IntegratorHard;
import etomica.potential.P1HardBoundary;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Osmotic pressure meter that calculates the difference in 
 */
public class MeterOsmoticPressure extends MeterPressureHard {
    
    private static final long serialVersionUID = 1L;
    private double virialSum;
    private double collisionRadius;
    private final P1HardBoundary[] leftBoundaries;
    private final P1HardBoundary[] rightBoundaries;
    
    public MeterOsmoticPressure(Space space, P1HardBoundary[] leftBoundaries, P1HardBoundary[] rightBoundaries) {
        super(space);
        this.leftBoundaries = leftBoundaries;
        this.rightBoundaries = rightBoundaries;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Records the number of simulation cycles performed by the integrator");
        return info;
    }
    
    /**
     * Sets the collision radius used to calculate accessible "area"
     * assuming that the relevant hard "boundaries" form right-angles with
     * other hard boundaries.  The given collision radius should be the 
     * collsion radius of the molecules with the adjacent boundaries.  If
     * the adjacent boundaries are periodic, this should be 0.
     */
    public void setCollisionRadius(double newRadius) {
        collisionRadius = newRadius;
    }
    
    /**
     * Returns the collision radius used to calculate the accessible "area".
     */
    public double getCollisionRadius() {
        return collisionRadius;
    }
    
    public void collisionAction(IntegratorHard.Agent agent) {
        for (int i=0; i<leftBoundaries.length; i++) {
            if (agent.collisionPotential == leftBoundaries[i]) {
                virialSum -= leftBoundaries[i].lastCollisionVirialTensor().component(0, 0);
            }
        }
        for (int i=0; i<rightBoundaries.length; i++) {
            if (agent.collisionPotential == rightBoundaries[i]) {
                virialSum += rightBoundaries[i].lastCollisionVirialTensor().component(0, 0);
            }
        }
    }

    public double getDataAsScalar() {
        double elapsedTime = timer.getDataAsScalar();
        double value = virialSum / elapsedTime;
        timer.reset();
        virialSum = 0;

        // calculate accessible "area"
        IVector dimensions = integratorHard.getPhase().getBoundary().getDimensions();
        double area = 1;
        for (int i=1; i<dimensions.getD(); i++) {
            area *= (dimensions.x(i)-2*collisionRadius);
        }
        return value / area;
    }
}