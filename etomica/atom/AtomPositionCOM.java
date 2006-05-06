package etomica.atom;

import java.io.Serializable;

import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.types.DataVector;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Length;

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

/*
 * History Created on Jan 26, 2005 by kofke
 */
public class AtomPositionCOM implements AtomPositionDefinition, Serializable {

    public AtomPositionCOM(Space space) {
        vectorSum = space.makeVector();
        center = space.makeVector();
        myAction = new MyAction(vectorSum);
        groupWrapper = new AtomGroupAction(myAction);
    }
    
    public Vector position(Atom atom) {
        vectorSum.E(0);
        myAction.massSum = 0;
        groupWrapper.actionPerformed(atom);
        center.Ea1Tv1(1.0 / myAction.massSum, vectorSum);
        return center;
    }

    private static class MyAction extends AtomActionAdapter {
        public MyAction(Vector sum) {
            vectorSum = sum;
        }
        public void actionPerformed(Atom a) {
            double mass = ((AtomTypeLeaf)a.type).getMass();
            vectorSum.PEa1Tv1(mass, ((AtomLeaf)a).coord.position());
            massSum += mass;
        }
        private Vector vectorSum;
        public double massSum;
    }
    
    private final Vector vectorSum;
    private final Vector center;
    private AtomGroupAction groupWrapper;
    private MyAction myAction;
}