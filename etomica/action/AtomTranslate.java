/*
 * History
 * Created on Nov 18, 2004 by kofke
 */
package etomica.action;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.Space;
import etomica.space.Vector;

public class AtomTranslate extends AtomActionAdapter {
    protected Vector displacement;
        
    public AtomTranslate(Space space) {
        super();
        displacement = space.makeVector();
    }
        
    public final void actionPerformed(Atom a) {((AtomLeaf)a).coord.position().PE(displacement);}
    public void actionPerformed(Atom a, Vector d) {((AtomLeaf)a).coord.position().PE(d);}
    public final void setDisplacement(Vector d) {displacement.E(d);}
}