/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica;


public class AtomTypeSphere extends AtomType {
    
    double diameter, radius;
    
    public AtomTypeSphere(AtomFactory creator) {
        this(creator, Default.ATOM_MASS, Default.ATOM_SIZE);
    }
    public AtomTypeSphere(AtomFactory creator, double m, double d) {
        super(creator, m);
        setDiameter(d);
    }
                
    public double diameter(Atom a) {return diameter;}
    public double radius(Atom a) {return radius;}
    
    /**
    * Sets diameter of this atom and updates radius accordingly.
    *
    * @param d   new value for diameter
    */
    public void setDiameter(double d) {diameter = d; radius = 0.5*d;}
}