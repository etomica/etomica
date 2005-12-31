package etomica.util;

import etomica.data.AccumulatorAverage;

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
            mass: amu
    */
    
    /**
     * Avogadro's number, 6.0221415e23.
     */
    //http://physics.nist.gov/cgi-bin/cuu/Value?na|search_for=abbr_in!
    public static final double AVOGADRO = 6.0221415e23;
        
    
    /**
     * The standard acceleration of gravity on Earth.
     */
//  acceleration of gravity (on Earth), in A/ps^2
    public static final double G = 9.8*1e10/1e24;  

    /**
     * Boltzmann's constant.
     */
//  Boltzmann's constant, converted from J/K to amu-A^2/ps^2 (where it equals 0.8314)
    public static final double BOLTZMANN_K = 1.380658e-23 * 1000 * AVOGADRO * 1e20 * 1e-24; 

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
    
    
    