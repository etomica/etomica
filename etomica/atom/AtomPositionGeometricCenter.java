package etomica.atom;

import etomica.action.AtomActionAdapter;
import etomica.action.AtomGroupAction;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * Calculates the geometric center over a set of atoms. The position of the
 * atom or child atoms are accumulated and used to compute their
 * center (unweighted by mass). Calculated center is obtained via the getPosition
 * method.
 * 
 * @author David Kofke
 */
public class AtomPositionGeometricCenter implements AtomPositionDefinition {

    public AtomPositionGeometricCenter(Simulation sim) {
        this(sim.space());
    }
    
    public AtomPositionGeometricCenter(Space space) {
        vectorSum = space.makeVector();
        center = space.makeVector();
        myAction = new MyAction(vectorSum);
        groupWrapper = new AtomGroupAction(myAction);
    }

    public Vector position(Atom atom) {
        vectorSum.E(0.0);
        myAction.nAtoms = 0;
        groupWrapper.actionPerformed(atom);
        center.Ea1Tv1(1.0 / myAction.nAtoms, vectorSum);
        return center;
    }

    private static class MyAction extends AtomActionAdapter {
        public MyAction(Vector v) {
            vectorSum = v;
            nAtoms = 0;
        }
        
        public void actionPerformed(Atom a) {
            vectorSum.PE(((AtomLeaf)a).getCoord().getPosition());
            nAtoms++;
        }
        
        private static final long serialVersionUID = 1L;
        private Vector vectorSum;
        public int nAtoms;
    }

    private final Vector center;
    private final Vector vectorSum;
    private AtomGroupAction groupWrapper;
    private MyAction myAction;
}