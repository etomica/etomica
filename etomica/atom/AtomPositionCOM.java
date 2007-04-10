package etomica.atom;

import java.io.Serializable;

import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Calculates the center of mass (COM) over a set of atoms. The mass and
 * position of all atoms passed to the actionPerformed are accumulated and used
 * to compute their center of mass. Calculated COM is obtained via the getData
 * method, which returns an array with the elements of the calculated COM
 * vector, or via the getCOM method, which returns a vector with the calculated
 * COM. Calculation is zeroed via the reset method.
 * <p>
 * A typical use of this class would have it passed to the allAtoms method of an
 * iterator, or wrapped in an AtomGroupAction to calculate the COM of the atoms
 * in an atom group.
 * 
 * @author David Kofke
 */

public class AtomPositionCOM implements AtomPositionDefinition, Serializable {

    public AtomPositionCOM(Space space) {
        vectorSum = space.makeVector();
        center = space.makeVector();
        myAction = new MyAction(vectorSum);
        groupWrapper = new AtomGroupAction(myAction);
    }
    
    public IVector position(Atom atom) {
        vectorSum.E(0);
        myAction.massSum = 0;
        groupWrapper.actionPerformed(atom);
        center.Ea1Tv1(1.0 / myAction.massSum, vectorSum);
        return center;
    }

    private static class MyAction extends AtomActionAdapter {
        public MyAction(IVector sum) {
            vectorSum = sum;
        }
        public void actionPerformed(Atom a) {
            double mass = ((AtomTypeLeaf)a.getType()).getMass();
            vectorSum.PEa1Tv1(mass, ((AtomLeaf)a).getPosition());
            massSum += mass;
        }
        private static final long serialVersionUID = 1L;
        private IVector vectorSum;
        public double massSum;
    }
    
    private static final long serialVersionUID = 1L;
    private final IVector vectorSum;
    private final IVector center;
    private AtomGroupAction groupWrapper;
    private MyAction myAction;
}