
package etomica.utility;
import java.io.*;

/**
 * Captures the functionality of the Java 1.2 class java.util.LinkedList
 * as used in the etomica api.
 *
 * Implements HashMap by wrapping a (Java 1.0) HashTable instance.  Consequently,
 * this implementation has only the features present in the Hashtable class, and
 * also is subject to the synchronization overhead inherent in Hashtable method calls.
 */

public class HashMap implements java.io.Serializable {

    private java.util.Hashtable hash;
    
    private static final Object NULL = new Object();
    
    public HashMap() {
        hash = new java.util.Hashtable();
    }
    public HashMap(int initialCapacity) {
        hash = new java.util.Hashtable(initialCapacity);
    }

    //used to facilitate null entries in hashmap
    private static final Object filterNull(Object obj) {
        return (obj == null) ? NULL : obj;
    }
    //used to facilitate null entries in hashmap
    private static final Object translate(Object obj) {
        return (obj == NULL) ? null : obj;
    }
    
    public void clear() {hash.clear();}
    
    public boolean containsKey(Object key) {
        return hash.containsKey(filterNull(key));
    }
    
    public boolean containsValue(Object value) {
        return hash.contains(filterNull(value));}
    
    public Object get(Object key) {return translate(hash.get(filterNull(key)));}
    
    public boolean isEmpty() {return hash.isEmpty();}
    
    public Object put(Object key, Object value) {return 
        translate(hash.put(filterNull(key), filterNull(value)));}
    
    public Object remove(Object key) {return 
        translate(hash.remove(filterNull(key)));}
    
    public int size() {return hash.size();}
    
    public Iterable keySet() {return new KeySet();}
    
    public Iterable values() {return new ValueSet();}
    
    
    public class KeySet implements Iterable {
        public Iterator iterator() {
            return new EnumerationWrapper(hash.keys());
        }
    }
    public class ValueSet implements Iterable {
        public Iterator iterator() {
            return new EnumerationWrapper(hash.elements());
        }
    }
    
    private class EnumerationWrapper implements Iterator {
        private final java.util.Enumeration enumeration;
        public EnumerationWrapper(java.util.Enumeration enum) {
            enumeration = enum;
        }
        public boolean hasNext() {return enumeration.hasMoreElements();}
        public Object next() {return translate(enumeration.nextElement());}
    }//end of EnumerationWrapper
    
        
}//end of HashMap