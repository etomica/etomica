package etomica.util;

import etomica.data.AccumulatorAverage;

/**
 * Collection of assorted physical constants.  All values
 * are in simulations units in which time is in picoseconds, length is in Angstroms,
 * and mass is in Daltons (amu).  Also defines several enumerated constants (or typed
 * constants).
 */
public interface Constants {
    
    public String VERSION = "01.11.20";    
    
    public static final double TWO_PI = 2.0*Math.PI;
    
    /*             
       Units for internal calculations (simulation units)
            time: ps
          length: Angstroms
            mass: amu
    */
    
    /**
     * Avogadro's number.
     */
    public static final double AVOGADRO = 6.0221367e23;
        
    
    /**
     * The standard acceleration of gravity on Earth.
     */
    public static final double G = 9.8*1e10/1e24;  //acceleration of gravity (on Earth), in A/ps^2

    /**
     * Boltzmann's constant.
     */
    public static final double BOLTZMANN_K = 1.380658e-23 * 1000 * AVOGADRO * 1e20 * 1e-24; //Boltzmann's constant, converted from J/K to amu-A^2/ps^2 (where it equals 0.8314)

    /**
     * TypedConstant classes are used to define a fixed set of specific values that can be taken by a field.
     * For example,  north/south/east/west.  Subclasses of this abstract class declare the general category
     * of the typed constant (e.g., Direction), and subclasses of the general class declare the set of allowable values.
     * Thus a field of type "Direction" may take take on only the values of the defined static final instances, namely
     * NORTH, SOUTH, EAST, WEST.  Further instances cannot be made because the constructor is private.
     * Fields may be compared against the static final values to query their value.  If the field is named "direction",
     * the construction is <code>if(direction == Constants.NORTH)</code>.
     * The constructor of a TypedConstant takes a String parameter that becomes the return value of the toString method.
     */
    public static abstract class TypedConstant implements java.io.Serializable {
        private final String label;
        protected TypedConstant(String s) {label = s;}
        public String toString() {return label;}
        public abstract TypedConstant[] choices();
    }
    
    /**
     * Typed constant for the directions TOP, BOTTOM, LEFT, RIGHT, FRONT, BACK.
     * Used to express the orientation of an object.
     */
     //maybe should rename this "Position", and use UP, DOWN, etc. for Direction
    public static class Direction extends TypedConstant {
        private Direction(String label) {super(label);}
        public static final Direction[] CHOICES = new Direction[] {
            new Direction("Top"),
            new Direction("Bottom"),
            new Direction("Left"),
            new Direction("Right"),
            new Direction("Front"),
            new Direction("Back")
        };
        public final TypedConstant[] choices() {return CHOICES;}
    }//end of Direction
    public static final Direction TOP = Direction.CHOICES[0];
    public static final Direction BOTTOM = Direction.CHOICES[1];
    public static final Direction LEFT = Direction.CHOICES[2];
    public static final Direction RIGHT = Direction.CHOICES[3];
    public static final Direction FRONT = Direction.CHOICES[4];
    public static final Direction BACK = Direction.CHOICES[5];
    
    /**
     * Typed constant for the compass directions NORTH, SOUTH, EAST, WEST.
     * Used to express the orientation of an object.
     */
    public static class CompassDirection extends TypedConstant {
        private CompassDirection(String label) {super(label);}
        public static final CompassDirection[] CHOICES = new CompassDirection[] {
            new CompassDirection("North"),
            new CompassDirection("South"),
            new CompassDirection("West"),
            new CompassDirection("East"),
        };
        public final TypedConstant[] choices() {return CHOICES;}
    }//end of CompassDirection
    public static final CompassDirection NORTH = CompassDirection.CHOICES[0];
    public static final CompassDirection SOUTH = CompassDirection.CHOICES[1];
    public static final CompassDirection WEST = CompassDirection.CHOICES[2];
    public static final CompassDirection EAST = CompassDirection.CHOICES[3];
    /**
     * Typed constant for specifying HORIZONTAL/VERTICAL alignment.
     */
    public static class Alignment extends TypedConstant {
        private Alignment(String label) {super(label);}
        public static final Alignment[] CHOICES = new Alignment[] {
            new Alignment("Horizontal (X Plane)"),
            new Alignment("Vertical (Y Plane)"),
            new Alignment("Width (Z Plane)")
        };
        public final TypedConstant[] choices() {return CHOICES;}
    }
    public static final Alignment HORIZONTAL = Alignment.CHOICES[0];
    public static final Alignment VERTICAL = Alignment.CHOICES[1];
    public static final Alignment WIDTH = Alignment.CHOICES[2];
    
    //MeterAbstract typed constants, repeated here to enable access by
    //implementing Constants interface
    public static final AccumulatorAverage.Type AVERAGE = AccumulatorAverage.AVERAGE;
    public static final AccumulatorAverage.Type ERROR = AccumulatorAverage.ERROR;
    public static final AccumulatorAverage.Type MOST_RECENT = AccumulatorAverage.MOST_RECENT;
    public static final AccumulatorAverage.Type MOST_RECENT_BLOCK = AccumulatorAverage.MOST_RECENT_BLOCK;
    public static final AccumulatorAverage.Type STANDARD_DEVIATION = AccumulatorAverage.STANDARD_DEVIATION;
    
}
    
    
    