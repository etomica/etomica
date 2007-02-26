package etomica.util;

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
    public static final double G_EARTH = 9.8*1e10/1e24;  

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
    
    /**
     * The speed of light, in simulation units.  Equal to 2997924.58 Angstroms/picosecond
     */
    public static final double LIGHT_SPEED = 299792458 * (1e10 * 1e-12);//convert m/s to A/ps
    
    /**
     * The gravitational constant, in simulation units.  
     * Equal to approximately 1.1e-31 A^3/ps^2/D. 
     */
    // 6.6742e-11 m^3 s^-2 kg^-1 (1e10 A/m)^3 (10-12 s/ps)^2 (1kg/1000g) (1g/Avo D)
    public static final double G = 6.6742e-11 * (1e30 * 1e-24 * 1e-3) /AVOGADRO;
    
    public static void main(String arg[]) {
        System.out.println("Avogadro's number: "+AVOGADRO);
        System.out.println("Boltzmann's constant: "+BOLTZMANN_K);
        System.out.println("C: "+LIGHT_SPEED);
        System.out.println("Planck's constant: "+PLANCK_H);
        System.out.println("Epsilon0: "+EPSILON_0);
        System.out.println("G: "+G);
        System.out.println("1.0/sqrt(4 Pi Epsilon0): "+1.0/Math.sqrt(4.*Math.PI*EPSILON_0));
        System.out.println("unit toSim: "+new etomica.units.systems.MKS().viscosity().toSim(1.0));
        System.out.println("unit toSim: "+etomica.units.Volt.UNIT.toSim(1.0));
        System.out.println("symbol: "+new LJ(1,1,1,false).viscosity().symbol());
    }
    /**
     * Enumerated type for the directions TOP, BOTTOM, LEFT, RIGHT, FRONT, BACK.
     * Used to express the orientation of an object.
     */
     //maybe should rename this "Position", and use UP, DOWN, etc. for Direction
    public static class Direction extends EnumeratedType {
        private Direction(String label) {super(label);}
        public static final Direction TOP = new Direction("Top");
        public static final Direction BOTTOM = new Direction("Bottom");
        public static final Direction LEFT = new Direction("Left");
        public static final Direction RIGHT = new Direction("Right");
        public static final Direction FRONT = new Direction("Front");
        public static final Direction BACK = new Direction("Back");
        public static Direction[] choices () {
            return new Direction[] {TOP,BOTTOM,LEFT,RIGHT,FRONT,BACK};
        }
        private static final long serialVersionUID = 1L;
    }
    
    /**
     * Enumerated type for the compass directions NORTH, SOUTH, EAST, WEST.
     * Used to express the orientation of an object.
     */
    public static class CompassDirection extends EnumeratedType {
        private CompassDirection(String label) {super(label);}
        public static final CompassDirection NORTH = new CompassDirection("North");
        public static final CompassDirection SOUTH = new CompassDirection("South");
        public static final CompassDirection WEST = new CompassDirection("West");
        public static final CompassDirection EAST = new CompassDirection("East");
        public static CompassDirection[] choicse() {
            return new CompassDirection[] {NORTH,SOUTH,WEST,EAST};
        }
        private static final long serialVersionUID = 1L;
   }
    /**
     * Enumerated type for specifying HORIZONTAL/VERTICAL alignment.
     */
    public static class Alignment extends EnumeratedType {
        private Alignment(String label) {super(label);}

        public static final Alignment HORIZONTAL = new Alignment("Horizontal (X Plane)");
        public static final Alignment VERTICAL = new Alignment("Vertical (Y Plane)");
        public static final Alignment WIDTH = new Alignment("Width (Z Plane)");

        public static Alignment[] choices() {
            return new Alignment[] {HORIZONTAL,VERTICAL,WIDTH};
        }
        private static final long serialVersionUID = 1L;
    }
    
}
    
    
