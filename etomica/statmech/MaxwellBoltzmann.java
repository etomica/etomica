package etomica.statmech;
import etomica.Simulation;

/**
 * Class with static methods for implementing features of the Maxwell-Boltzmann velocity distribution
 */
public final class MaxwellBoltzmann {
        
    /**
     * Private constructor to prevent instantiation
     */
    private MaxwellBoltzmann() {}
    
    /**
     * Returns a random velocity component for a particle of mass m at temperature T
     */
    public static double randomMomentumComponent(double T, double m) {
        return Simulation.random.nextGaussian()*Math.sqrt(m*T);
    }
}