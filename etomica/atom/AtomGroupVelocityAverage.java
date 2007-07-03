package etomica.atom;

import java.io.Serializable;

import etomica.action.AtomAction;
import etomica.action.AtomGroupAction;
import etomica.space.IVector;
import etomica.space.Space;

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
    public IVector getVelocityAverage(IAtom a) {
        vectorSum.E(0.0);
        myAction.massSum = 0;
        groupWrapper.actionPerformed(a);
        vectorSum.TE(1.0 / myAction.massSum);
        return vectorSum;
    }

    private static class MyAction implements AtomAction {
        public MyAction(IVector sum) {
            vectorSum = sum;
        }
        public void actionPerformed(IAtom a) {
            vectorSum.PE(((IAtomKinetic)a).getVelocity());
            massSum += ((AtomTypeLeaf)a.getType()).getMass();
        }
        private static final long serialVersionUID = 1L;
        private final IVector vectorSum;
        public double massSum;
    }

    private static final long serialVersionUID = 1L;
    private final IVector vectorSum;
    private MyAction myAction;
    private AtomGroupAction groupWrapper;
}