package etomica.utility;
import java.util.HashMap;

/**
 * Implementation of a hash map in which entries are keyed to a pair of objects,
 * instead of to a single object.
 */
public final class HashMap2 implements java.io.Serializable {
    
    public HashMap table;
    int defaultCapacity;
    
    public HashMap2() {
        this(100);
    }
    public HashMap2(int n) {
        defaultCapacity = n;
        table = new HashMap(n);
    }
    
    /**
     * Adds the given value to the map, associating it with the two keys.
     * Order of the keys is not relevant.  Value can be obtained by
     * call to either get(key1, key2) or get(key2, key1).
     */
    public Object put(Object key1, Object key2, Object value) {
        putOrdered(key2, key1, value);
        return putOrdered(key1, key2, value);
    }
    /**
     * Adds the given value to the map, associating it with the two keys.
     * Order of the keys is relevant.  Value can be obtained only
     * by call to get(key1, key2); call to get(key2, key1) returns different value.
     */
    public Object putOrdered(Object key1, Object key2, Object value) {
        HashMap subTable = get(key1);
        if(subTable == null) {
            subTable = new HashMap(defaultCapacity);
            table.put(key1,subTable);
        }
        return subTable.put(key2, value);
    }
    
    /**
     * Returns the objects associated with the two keys.
     */
    public Object get(Object key1, Object key2) {
        HashMap subTable = get(key1);
        return (subTable == null) ? null : subTable.get(key2);
    }
    
    /**
     * Returns a hash map for all values in which argument is one of the keys.
     */
     public HashMap get(Object key1) {
        return (HashMap)table.get(key1);
     }
    
    /**
     * main method to test and demonstrate this class
     */
    public static void main(String[] args) {
        
        HashMap2 table = new HashMap2();
        table.put("A","A","AA");
        table.put("A","a","Aa");
        table.put("B","a","Ba");
        table.put("B","A","BA");
        System.out.println(table.get("A","A"));//AA
        System.out.println(table.get("A","a"));//Aa
        System.out.println(table.get("B","a"));//Ba
        System.out.println(table.get("B","A"));//BA
        System.out.println(table.get("C","a"));//null
        System.out.println(table.get("A","C"));//null
    }
}