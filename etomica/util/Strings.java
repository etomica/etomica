package etomica.util;

/**
 * Non-instantiable class providing a few static methods for String manipulation.
 */
public class Strings {
    
    /**
     * Private construction to ensure non-instantiation.
     */
    private Strings() {}
    
    /**
     * Capitalizes the first letter of the given string.
     * Can accept zero-length or null argument, and simply returns the same.
     */
    public static String capitalize(String s) {
        if(s == null || s.length() == 0) return s;
        return s.substring(0,1).toUpperCase() + s.substring(1);
    }
    
    /**
     * Converts the first letter of the given string to lower case.
     * Can accept zero-length or null argument, and simply returns the same.
     */
    public static String decapitalize(String s) {
        if(s == null || s.length() == 0) return s;
        return s.substring(0,1).toLowerCase() + s.substring(1);
    }
    
    /**
     * Takes an array of objects and converts them to an array of strings,
     * with each string given by the object's toString method.
     */
    public static String[] toStringArray(Object[] obj) {
        String[] s = new String[obj.length];
        for(int i=0; i<obj.length; i++) {s[i] = obj[i].toString();}
        return s;
    }
    
    /**
     * Takes the elements of two string arrays and collects them
     * into a new array containing both sets of strings.
     */
    public static String[] arrayCollect(String[] s1, String[] s2) {
        String[] s = new String[s1.length + s2.length];
        for(int i=0; i<s1.length; i++) {s[i] = s1[i];}
        for(int i=0; i<s2.length; i++) {s[s1.length+i] = s2[i];}
        return s;
    }
    
    /**
     * Returns a String that formats the given value as an exponent.  If
     * the given value is equal to 1.0, returns an empty string.
     */
    //this method may be expanded to handle other types of formatting, such as HTML
    public static String exponent(double v) {
        if(v == 1) {
            return "";
        } else if ((int)v == v) {
            return "^"+(int)v;
        } else {
            return "^"+v;
        }
    }
    
    /**
     * Returns a String that formats the given String as an exponent.  If
     * the given value is equal to 1.0, returns an empty string.
     */
    public static String exponent(String v) {
        if(v.equals("")) {
            return v;
        } else {
            return "^"+v;
        }
    }
    
    public static void main(String[] args) {
        System.out.println("First argument capitalized: "+capitalize(args[0]));
        System.out.println("Second argument decapitalized: "+decapitalize(args[1]));
    }
}
    
    
