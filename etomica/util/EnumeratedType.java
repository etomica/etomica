package etomica.util;




/*
 * History
 * Created on Aug 24, 2005 by kofke
 */
/**
 * EnumeratedType classes are used to define a fixed set of specific values that can be taken by a field.
 * For example,  north/south/east/west.  Subclasses of this abstract class declare the general category
 * of the typed constant (e.g., Direction), and subclasses of the general class declare the set of allowable values.
 * Thus a field of type "Direction" may take take on only the values of the defined static final instances, namely
 * NORTH, SOUTH, EAST, WEST.  Further instances cannot be made because the constructor is private.
 * Fields may be compared against the static final values to query their value.  If the field is named "direction",
 * the construction is <code>if(direction == Constants.NORTH)</code>.
 * The constructor of a EnumeratedType takes a String parameter that becomes the return value of the toString method.
 */
public abstract class EnumeratedType implements java.io.Serializable {
    private final String label;
    protected EnumeratedType(String s) {label = s;}
    public String toString() {return label;}
}