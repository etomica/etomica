package etomica.statmech;
import etomica.Simulation;

/**
 * Class with static methods for implementing features of the Maxwell-Boltzmann velocity distribution
 */
 
 /* History
  * 09/08/02 (DAK) implemented as function
  */
  
public final class MaxwellBoltzmann {
    
    //private to prevent instantiation
    private MaxwellBoltzmann() {}
        
    /**
     * Returns a random velocity component for a particle of mass m at temperature T
     */
    public static double randomMomentumComponent(double T, double m) {
        return Simulation.random.nextGaussian()*Math.sqrt(m*T);
    }
    
    /**
     * Function giving the Maxwell-Boltzmann distribution of speeds.
     */
    public static class Distribution implements etomica.utility.Function {
        private double temperature;
        private double mass;
        private double a, c;
        public Distribution() {
            this(etomica.Default.TEMPERATURE, etomica.Default.ATOM_MASS);
        }
        public Distribution(double temperature, double mass) {
            this.temperature = temperature;
            this.mass = mass;
            update();
        }
        public double f(double v) {
            double v2 = v*v;
            return c*v2*Math.exp(-a*v2);
        }
        //not well checked.
        public double dfdx(double v) {
            double v2 = v*v;
            return c*Math.exp(-a*v2)*2*v*(1 - a*v2);
        }
        public double inverse(double x) {
            throw new RuntimeException("MaxwellBoltxmann.Distribution.inverse not defined");
        }
        public void setTemperature(double t) {
            temperature = t;
            update();
            
        }
        private void update() {
            c = 4.0*Math.PI*Math.pow(0.5*mass/Math.PI/temperature,1.5);
            a = 0.5*mass/temperature;
        }
        public double getTemperature() {return temperature;}
        public etomica.units.Dimension getTemperatureDimension() {return etomica.units.Dimension.TEMPERATURE;}

        public void setMass(double m) {
            mass = m;
            update();
        }
        public double getMass() {return mass;}
        public etomica.units.Dimension getMassDimension() {return etomica.units.Dimension.MASS;}
    }//end of Distribution
        
    
    
}