package etomica.atom;

import java.io.Serializable;

import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Calculates the mass average velocity over a set of atoms. The velocity
 * and mass of the atom passed to getVelocityAverage all its child atoms 
 * are used to compute the mass-average velocity (total momentum divided by 
 * the total mass).
 * 
 * @author David Kofke
 */
public class AtomGroupVelocityAverage implements Serializable {

    public AtomGroupVelocityAverage(Space space) {
        vectorSum = space.makeVector();
        myAction = new MyAction(vectorSum);
        groupWrapper = new AtomGroupAction(myAction);
    }
    
    /**
     * Returns the mass-average velocity of the given Atom and 
     * all its children.
     */
    public Vector getVelocityAverage(Atom a) {
        vectorSum.E(0.0);
        myAction.massSum = 0;
        groupWrapper.actionPerformed(a);
        vectorSum.TE(1.0 / myAction.massSum);
        return vectorSum;
    }

    private static class MyAction extends AtomActionAdapter {
        public MyAction(Vector sum) {
            vectorSum = sum;
        }
        public void actionPerformed(Atom a) {
            vectorSum.PE(((ICoordinateKinetic)((AtomLeaf)a).coord).velocity());
            massSum += ((AtomTypeLeaf)a.type).getMass();
        }
        private final Vector vectorSum;
        public double massSum;
    }

    private final Vector vectorSum;
    private MyAction myAction;
    private AtomGroupAction groupWrapper;
}