package simulate.utility;

public class StringUtility {
    
    public static String capitalize(String s) {
        return s.substring(0,1).toUpperCase() + s.substring(1);
    }
    
    public static String decapitalize(String s) {
        return s.substring(0,1).toLowerCase() + s.substring(1);
    }
    
    /**
     * Takes and array of objects and converts them to an array of strings,
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
    
    
    public static void main(String[] args) {
        System.out.println("First argument capitalized: "+capitalize(args[0]));
        System.out.println("Second argument decapitalized: "+decapitalize(args[1]));
    }
}
    
    