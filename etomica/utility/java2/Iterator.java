
package etomica.utility.java2;

/* History
 *  
 * 12/31/02 (DAK) appended java2 to package name
 * 07/28/04 (DAK) added reset()
 */
public interface Iterator {
    boolean hasNext();
    Object next();
    void reset();
}
