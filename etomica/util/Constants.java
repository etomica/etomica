package etomica.util;

import etomica.data.AccumulatorAverage;
import etomica.units.Joule;
import etomica.units.systems.LJ;

/**
 * Collection of assorted physical constants.  All values
 * are in simulations units in which time is in picoseconds, length is in Angstroms,
 * and mass is in Daltons (amu).  Also defines several enumerated-type constants.
 */
public final class Constants {

    //private constructor to prevent instantiation
    private Constants() {
    }
    
    public static final double TWO_PI = 2.0*Math.PI;
    
    /*             
       Units for internal calculations (simulation units)
            time: ps
          length: Angstroms
            mass: Daltons, or atomic mass units (amu)
    */
    
    /**
     * Avogadro's number, 6.0221415e23 molecules/mole.
     */
    // reference:  http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=abbr_in!
    public static final double AVOGADRO = 6.0221415e23;
        
    
    /**
     * The standard acceleration of gravity on Earth.
     */
    //  acceleration of gravity (on Earth), in A/ps^2
    public static final double G = 9.8*1e10/1e24;  

    /**
     * Boltzmann's constant, in (sim units)/kelvin.  Specifically,
     * equal to 0.8314517839107 (daltons)(angstroms^2)(ps^-2)(kelvin^-1).
     */
    //  Boltzmann's constant, converted from J/K to D-A^2/ps^2/K (where it equals 0.8314)
    // (1.38e-23 kg-m^2/s^2/K/molecule)(1000 g/kg)(N_avo D/g)(10^10 A/m)^2 (10^-12 s/ps)^2  
    public static final double BOLTZMANN_K = 1.380658e-23 * 1000 * AVOGADRO * 1e20 * 1e-24; 

    /**
     * Planck's constant, in simulation units.  Specifically, equal to approximately 39.903127 D-A^2/ps.
     */
    //convert from J-s to simulation units
    public static final double PLANCK_H = Joule.UNIT.toSim(6.6260693e-34) * 1e12;
    
    /**
     * The permittivity of a vacuum, in (sim units) (electron charge)^2. Specifically,
     * equal to  5.7276575390366254E-7 (ps^2)(daltons^-1)(angstroms^-3)(electrons^2).
     */
    //epsilon0, converted from C^2/(N-m^2) to e^2 ps^2/(D-A^3)
    // (8.854e-12 C^2 s^2/(kg-m^3)) (1/1.60217653e-12 e/C)^2 (10^12 ps/s)^2 (10^-3 kg/g) (1/Avo g/D) (10^-10 m/A)^3
    public static final double EPSILON_0 = 8.8541878176e-12 * (1.0/1.60217653e-19/1.60217653e-19)
            * 1e24 * 1e-3 / AVOGADRO * 1e-30;
    
    public static final double LIGHT_SPEED = 299792458 * 1e10 * 1e-12;//convert m/s to A/ps
    
    public static void main(String arg[]) {
        System.out.println("Avogadro's number: "+AVOGADRO);
        System.out.println("Boltzmann's constant: "+BOLTZMANN_K);
        System.out.println("Planck's constant: "+PLANCK_H);
        System.out.println("Epsilon0: "+EPSILON_0);
        System.out.println("1.0/sqrt(4 Pi Epsilon0): "+1.0/Math.sqrt(4.*Math.PI*EPSILON_0));
        System.out.println("unit toSim: "+etomica.units.systems.MKS.SYSTEM.viscosity().toSim(1.0));
        System.out.println("unit toSim: "+etomica.units.Poise.UNIT.toSim(1.0));
        System.out.println("symbol: "+new LJ(1,1,1,false).viscosity().symbol());
    }
    /**
     * Enumerated type for the directions TOP, BOTTOM, LEFT, RIGHT, FRONT, BACK.
     * Used to express the orientation of an object.
     */
     //maybe should rename this "Position", and use UP, DOWN, etc. for Direction
    public static class Direction extends EnumeratedType {
        private Direction(String label) {super(label);}
        public static final Direction[] CHOICES = new Direction[] {
            new Direction("Top"),
            new Direction("Bottom"),
            new Direction("Left"),
            new Direction("Right"),
            new Direction("Front"),
            new Direction("Back")
        };
        public final EnumeratedType[] choices() {return CHOICES;}
    }//end of Direction
    public static final Direction TOP = Direction.CHOICES[0];
    public static final Direction BOTTOM = Direction.CHOICES[1];
    public static final Direction LEFT = Direction.CHOICES[2];
    public static final Direction RIGHT = Direction.CHOICES[3];
    public static final Direction FRONT = Direction.CHOICES[4];
    public static final Direction BACK = Direction.CHOICES[5];
    
    /**
     * Enumerated type for the compass directions NORTH, SOUTH, EAST, WEST.
     * Used to express the orientation of an object.
     */
    public static class CompassDirection extends EnumeratedType {
        private CompassDirection(String label) {super(label);}
        public static final CompassDirection[] CHOICES = new CompassDirection[] {
            new CompassDirection("North"),
            new CompassDirection("South"),
            new CompassDirection("West"),
            new CompassDirection("East"),
        };
        public final EnumeratedType[] choices() {return CHOICES;}
    }//end of CompassDirection
    public static final CompassDirection NORTH = CompassDirection.CHOICES[0];
    public static final CompassDirection SOUTH = CompassDirection.CHOICES[1];
    public static final CompassDirection WEST = CompassDirection.CHOICES[2];
    public static final CompassDirection EAST = CompassDirection.CHOICES[3];
    /**
     * Enumerated type for specifying HORIZONTAL/VERTICAL alignment.
     */
    public static class Alignment extends EnumeratedType {
        private Alignment(String label) {super(label);}
        public static final Alignment[] CHOICES = new Alignment[] {
            new Alignment("Horizontal (X Plane)"),
            new Alignment("Vertical (Y Plane)"),
            new Alignment("Width (Z Plane)")
        };
        public final EnumeratedType[] choices() {return CHOICES;}
    }
    public static final Alignment HORIZONTAL = Alignment.CHOICES[0];
    public static final Alignment VERTICAL = Alignment.CHOICES[1];
    public static final Alignment WIDTH = Alignment.CHOICES[2];
    
    public static final AccumulatorAverage.Type AVERAGE = AccumulatorAverage.AVERAGE;
    public static final AccumulatorAverage.Type ERROR = AccumulatorAverage.ERROR;
    public static final AccumulatorAverage.Type MOST_RECENT = AccumulatorAverage.MOST_RECENT;
    public static final AccumulatorAverage.Type MOST_RECENT_BLOCK = AccumulatorAverage.MOST_RECENT_BLOCK;
    public static final AccumulatorAverage.Type STANDARD_DEVIATION = AccumulatorAverage.STANDARD_DEVIATION;
    
}
    
    
    