/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.atom;

import etomica.Atom;
import etomica.AtomIndexManager;
import etomica.AtomType;
import etomica.Default;


public class AtomTypeSphere extends AtomType {
    
    double diameter, radius;
    
    public AtomTypeSphere(AtomIndexManager indexManager) {
        this(indexManager, Default.ATOM_MASS, Default.ATOM_SIZE);
    }
    public AtomTypeSphere(AtomIndexManager indexManager, double m, double d) {
        super(indexManager, m);
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