package etomica.normalmode;

import java.io.Serializable;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.space.IVector;

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

    public CoordinateDefinitionLeaf(Phase phase, Primitive primitive) {
        this(phase, primitive, new BasisMonatomic(phase.getSpace()));
    }
    
    public CoordinateDefinitionLeaf(Phase phase, Primitive primitive, Basis basis) {
        super(phase, phase.getSpace().D()*basis.getScaledCoordinates().length, primitive, basis);
        workVector = phase.getSpace().makeVector();
        u = new double[coordinateDim];
    }

    /**
     * Assigns the given array u to be the current position of the atom minus its lattice position
     */
    public double[] calcU(AtomSet atoms) {
        int j = 0;
        for (int i=0; i<atoms.getAtomCount(); i++) {
            IAtomPositioned a = (IAtomPositioned)atoms.getAtom(i);
            IVector pos = a.getPosition();
            IVector site = getLatticePosition(a);
            workVector.Ev1Mv2(pos, site);
            for (int k=0; k<workVector.getD(); k++) {
                u[j+k] = workVector.x(k);
            }
            j += workVector.getD();
        }
        return u;
    }

    public void initNominalU(AtomSet molecules) {
        //nothing to do -- lattice site is all information needed for u
    }

    /**
     * Sets the position of the atom to be its lattice position plus the offset u
     */
    public void setToU(AtomSet atoms, double[] newU) {
        int j = 0;
        for (int i=0; i<atoms.getAtomCount(); i++) {
            IAtomPositioned a = (IAtomPositioned)atoms.getAtom(i);
            IVector pos = a.getPosition();
            for (int k=0; k<workVector.getD(); k++) {
                pos.setX(k, newU[j+k]);
            }
            j += workVector.getD();
            pos.PE(getLatticePosition(a));
        }
    }

    protected final IVector workVector;
    protected final double[] u;
    private static final long serialVersionUID = 1L;
}
