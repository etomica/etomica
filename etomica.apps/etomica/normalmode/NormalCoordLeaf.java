package etomica.normalmode;

import java.io.Serializable;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * NormalCoordWrapper implementation for leaf atoms.  The class simply takes
 * the u values to be real space displacements from the nominal positions.
 * @author Andrew Schultz
 */
public class NormalCoordLeaf implements NormalCoordMapper, Serializable {

    public NormalCoordLeaf(Space space) {
        this.space = space;
    }

    public void setNumAtoms(int numAtoms) {
        nominalU = new double[numAtoms][getNormalDim()];
    }
    
    public int getNormalDim() {
        return space.D();
    }

    public void calcU(Atom atom, int atomCount, double[] u) {
        Vector pos = ((AtomLeaf)atom).coord.position();
        for (int i=0; i<pos.D(); i++) {
            u[i] = pos.x(i) - nominalU[atomCount][i];
        }
    }

    public void initNominalU(Atom atom, int atomCount) {
        Vector pos = ((AtomLeaf)atom).coord.position();
        for (int i=0; i<pos.D(); i++) {
            nominalU[atomCount][i] = pos.x(i);
        }
    }
    
    public void setToU(Atom atom, int atomCount, double[] u) {
        ((AtomLeaf)atom).coord.position().E(u);
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected double[][] nominalU;
}
