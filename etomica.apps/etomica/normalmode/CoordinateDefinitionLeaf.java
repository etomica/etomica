package etomica.normalmode;

import java.io.Serializable;

import etomica.atom.Atom;
import etomica.atom.AtomLeaf;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * CoordinateDefinition implementation for monatomic molecules that are simply
 * leaf atoms. The class simply takes the u values to be real space
 * displacements from the nominal positions.
 * 
 * @author Andrew Schultz
 */
public class CoordinateDefinitionLeaf extends CoordinateDefinition implements
        Serializable {

    public CoordinateDefinitionLeaf(Space space) {
        super(space.D());
    }

    public void calcU(Atom atom, int index, double[] u) {
        IVector pos = ((AtomLeaf) atom).getCoord().getPosition();
        for (int i = 0; i < pos.getD(); i++) {
            u[i] = pos.x(i) - nominalU[index][i];
        }
    }

    public void initNominalU(Atom atom, int index) {
        IVector pos = ((AtomLeaf) atom).getCoord().getPosition();
        for (int i = 0; i < pos.getD(); i++) {
            nominalU[index][i] = pos.x(i);
        }
    }

    public void setToU(Atom atom, int index, double[] u) {
        IVector pos = ((AtomLeaf) atom).getCoord().getPosition();
        for (int i = 0; i < pos.getD(); i++) {
            pos.setX(i, nominalU[index][i] + u[i]);
        }
    }

    public void setNumAtoms(int numAtoms) {
        nominalU = new double[numAtoms][getCoordinateDim()];
    }

    private static final long serialVersionUID = 1L;
}
