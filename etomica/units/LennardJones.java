package simulate.units;
import simulate.Default;

/**
 * Lennard-Jones units.
 * All quantities are made dimensionless with respect to a characteristic
 * size (sigma), energy (epsilon) and mass.
 * Values of sigma, epsilon and mass may be set to any value using the static accessor
 * methods (set/get).  Changes are immediately propagated to the sim->LJ unit conversion factors.
 */
 
public class LennardJones implements java.io.Serializable {
    
    private static double sigma = Default.ATOM_SIZE;
    private static double epsilon = Default.POTENTIAL_WELL;
    private static double mass = Default.ATOM_MASS;
        
    public static double getSigma() {return sigma;}
    public static void setSigma(double s) {sigma = s; update();}
    public static double getEpsilon() {return epsilon;}
    public static void setEpsilon(double e) {epsilon = e; update();}
    public static double getMass() {return mass;}
    public static void setMass(double m) {mass = m; update();}
        
    private static void update() {
        Mass.UNIT.updateConversions();
        Length.UNIT.updateConversions();
        Time.UNIT.updateConversions();
        Charge.UNIT.updateConversions();
        Dipole.UNIT.updateConversions();
        Energy.UNIT.updateConversions();
        Temperature.UNIT.updateConversions();
        Pressure.UNIT.updateConversions();
        Pressure2D.UNIT.updateConversions();
        Volume.UNIT.updateConversions();
        Volume2D.UNIT.updateConversions();
    }
    
    //  \u03B5 is unicode for epsilon
    //  \u03C3 is unicode for sigma
    //  \u00BD is unicode for "1/2"
    
    public static final class Mass extends BaseUnit.Mass {
        public static final LennardJones.Mass UNIT = new LennardJones.Mass();
        public Mass() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "m";   
            updateConversions();
        }
        public void updateConversions() {
            to = mass; //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Length extends BaseUnit.Length {
        public static final LennardJones.Length UNIT = new LennardJones.Length();
        public Length() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03C3";   
            updateConversions();
        }
        public void updateConversions() {
            to = sigma; //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Time extends BaseUnit.Time {
        public static final LennardJones.Time UNIT = new LennardJones.Time();
        public Time() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "(\u03B5m\u03C3^2)^\u00BD";   
            updateConversions();
        }
        public void updateConversions() {
            to = 1.0/Math.sqrt(epsilon * mass * sigma*sigma); //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Charge extends BaseUnit.Charge {
        public static final LennardJones.Charge UNIT = new LennardJones.Charge();
        public Charge() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "(\u03B5\u03C3)^(-\u00BD)";   
            updateConversions();
        }
        public void updateConversions() {
            to = Math.sqrt(epsilon * sigma); //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Dipole extends BaseUnit.Dipole {
        public static final LennardJones.Dipole UNIT = new LennardJones.Dipole();
        public Dipole() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "(\u03B5\u03C3^3)^(-\u00BD)";   
            updateConversions();
        }
        public void updateConversions() {
            to = Math.sqrt(epsilon * sigma*sigma*sigma); //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Energy extends BaseUnit.Energy {
        public static final LennardJones.Energy UNIT = new LennardJones.Energy();
        public Energy() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03B5";   
            updateConversions();
        }
        public void updateConversions() {
            to = epsilon; //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }


    public static final class Temperature extends BaseUnit.Energy {
        public static final LennardJones.Temperature UNIT = new LennardJones.Temperature();
        public Temperature() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03B5";   
            updateConversions();
        }
        public void updateConversions() {
            to = epsilon; //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }
    public static final class Pressure extends BaseUnit.Pressure {
        public static final LennardJones.Pressure UNIT = new LennardJones.Pressure();
        public Pressure() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03B5/\u03C3^3";   
            updateConversions();
        }
        public void updateConversions() {
            to = epsilon/(sigma*sigma*sigma); //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Pressure2D extends BaseUnit.Pressure2D {
        public static final LennardJones.Pressure2D UNIT = new LennardJones.Pressure2D();
        public Pressure2D() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03B5/\u03C3^2";   
            updateConversions();
        }
        public void updateConversions() {
            to = epsilon/(sigma*sigma); //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Volume extends BaseUnit.Volume {
        public static final LennardJones.Volume UNIT = new LennardJones.Volume();
        public Volume() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03C3^3";   
            updateConversions();
        }
        public void updateConversions() {
            to = sigma*sigma*sigma; //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }

    public static final class Volume2D extends BaseUnit.Volume2D {
        public static final LennardJones.Volume2D UNIT = new LennardJones.Volume2D();
        public Volume2D() {
            prefixAllowed = false;
            name = "LJ units";
            symbol = "\u03C3^2";   
            updateConversions();
        }
        public void updateConversions() {
            to = sigma*sigma; //conversion from LJ unit to simulation unit
            from = 1.0/to;
        }
    }
}