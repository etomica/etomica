package etomica.normalmode;

import java.io.Serializable;

import etomica.atom.AtomLeaf;
import etomica.atom.IAtom;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * CoordinateDefinition implementation for monatomic molecules that are simply
 * leaf atoms. The class simply takes the u values to be real space
 * displacements from the nominal positions.
 * 
 * @author Andrew Schultz
 */

//when dealing with heterogeneous molecular systems we may need to introduce the mass in the definition

public class CoordinateDefinitionLeaf extends CoordinateDefinition implements
        Serializable {

    public CoordinateDefinitionLeaf(Space space) {
        super(space.D());
        workVector = space.makeVector();
    }

    /**
     * Assigns the given array u to be the current position of the atom minus its lattice position
     */
    public void calcU(IAtom[] atom, double[] u) {
        IVector pos = ((AtomLeaf) atom[0]).getPosition();
        IVector site = getLatticePosition(atom[0]);
        workVector.Ev1Mv2(pos, site);
        workVector.assignTo(u);
    }

    public void initNominalU(IAtom[] atom) {
        //nothing to do -- lattice site is all information needed for u
    }

    /**
     * Sets the position of the atom to be its lattice position plus the offset u
     */
    public void setToU(IAtom[] atom, double[] u) {
        workVector.E(u);
        IVector site = getLatticePosition(atom[0]);
        ((AtomLeaf) atom[0]).getPosition().Ev1Pv2(site, workVector);
    }

    public void setNumAtoms(int numAtoms) {
        nominalU = new double[numAtoms][getCoordinateDim()];
    }

    private final IVector workVector;
    private static final long serialVersionUID = 1L;
}
