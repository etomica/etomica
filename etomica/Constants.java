package etomica;
import java.awt.Color;
import java.util.Random;

public class Constants extends Object {
    
    public String getVersion() {return "01.03.10";}    
    private static final Random random = new Random();
    
    private Constants() {}   // can't instantiate class
    
    public static final double TWO_PI = 2.0*Math.PI;
    
    /*             
       Units for internal calculations (simulation units)
            time: ps
          length: Angstroms
            mass: amu
    */
    
    public static final double AVOGADRO = 6.0221367e23;
        
    
    public static final double G = 9.8*1e10/1e24;  //acceleration of gravity (on Earth), in A/ps^2
    public static final double BOLTZMANN_K = 1.380658e-23 * 1000 * AVOGADRO * 1e20 * 1e-24; //Boltzmann's constant, converted from J/K to amu-A^2/ps^2 (where it equals 0.8314)

  /* Colors adopted in the web textbook on molecular simulation */
    public static final Color KHAKI = new Color(153,153,102);
    public static final Color DARK_KHAKI = new Color(102,102,051);
    public static final Color BRIGHT_RED = new Color(153,000,000);
    public static final Color DARK_RED = new Color(102,000,000);
    public static final Color BLUSH = new Color(153,102,102);
    public static final Color TAN = new Color(204,204,153);
    public static final Color RandomColor() {return new Color(random.nextFloat(),random.nextFloat(),random.nextFloat());}
    public static final Color randomColor() {return new Color(random.nextFloat(),random.nextFloat(),random.nextFloat());}
  
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
     * Typed constant for the compass directions NORTH, SOUTH, EAST, WEST.
     * Used to express the orientation of an object.
     */
    public static class Direction extends TypedConstant {
        private Direction(String label) {super(label);}
        public static final Direction[] CHOICES = new Direction[] {
            new Direction("North"),
            new Direction("East"),
            new Direction("South"),
            new Direction("West")
        };
        public final TypedConstant[] choices() {return CHOICES;}
    }//end of Direction
    public static final Direction NORTH = Direction.CHOICES[0];
    public static final Direction EAST = Direction.CHOICES[1];
    public static final Direction SOUTH = Direction.CHOICES[2];
    public static final Direction WEST = Direction.CHOICES[3];
    
    /**
     * Typed constant for specifying HORIZONTAL/VERTICAL alignment.
     */
    public static class Alignment extends TypedConstant {
        private Alignment(String label) {super(label);}
        public static final Alignment[] CHOICES = new Alignment[] {
            new Alignment("Horizontal"),
            new Alignment("Vertical")};
        public final TypedConstant[] choices() {return CHOICES;}
    }
    public static final Alignment HORIZONTAL = Alignment.CHOICES[0];
    public static final Alignment VERTICAL = Alignment.CHOICES[1];
}
    
    
    