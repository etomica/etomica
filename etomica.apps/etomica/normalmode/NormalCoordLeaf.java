package etomica.normalmode;

import java.io.Serializable;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.IVector;
import etomica.space.Space;

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
        IVector pos = ((AtomLeaf)atom).getCoord().getPosition();
        for (int i=0; i<pos.getD(); i++) {
            u[i] = pos.x(i) - nominalU[atomCount][i];
        }
    }

    public void initNominalU(Atom atom, int atomCount) {
        IVector pos = ((AtomLeaf)atom).getCoord().getPosition();
        for (int i=0; i<pos.getD(); i++) {
            nominalU[atomCount][i] = pos.x(i);
        }
    }
    
    public void setToU(Atom atom, int atomCount, double[] u) {
        IVector pos = ((AtomLeaf)atom).getCoord().getPosition();
        for (int i=0; i<pos.getD(); i++) {
            pos.setX(i,nominalU[atomCount][i]+u[i]);
        }
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected double[][] nominalU;
}
