package etomica.normalmode;

import java.io.Serializable;

import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.box.Box;
import etomica.space.IVector;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz
 */
public class CoordinateDefinitionMolecule extends CoordinateDefinition
        implements Serializable {

    public CoordinateDefinitionMolecule(Box box, Primitive primitive, int orientationDim) {
        this(box, primitive, orientationDim, new BasisMonatomic(box.getSpace()));
    }
    
    public CoordinateDefinitionMolecule(Box box, Primitive primitive, int orientationDim, Basis basis) {
        super(box, (box.getSpace().D() + orientationDim)*basis.getScaledCoordinates().length, primitive, basis);
        work1 = box.getSpace().makeVector();
        u = new double[coordinateDim];
    }

    public double[] calcU(AtomSet molecules) {
        // calculates components of U related to the the center of mass of the
        // molecules
        // subclass is responsible for setting orientation or intramolecular
        // degrees of freedom
        int j = 0;
        for (int i=0; i<molecules.getAtomCount(); i++) {
            IAtom molecule = molecules.getAtom(i);
            IVector pos = molecule.getType().getPositionDefinition().position(molecule);
            IVector site = getLatticePosition(molecule);
            work1.Ev1Mv2(pos, site);
            for (int k = 0; k < pos.getD(); k++) {
                u[j+k] = work1.x(k);
            }
            j += coordinateDim/molecules.getAtomCount();

        }
        return u;
    }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(AtomSet molecules) {
    }
    
    public int getULength(){
        return u.length;
    }

    public void setToU(AtomSet molecules, double[] newU) {
        // sets the center of mass of the molecules to that specified by newU
        // subclass is responsible for setting orientation or intramolecular
        // degrees of freedom
        int j = 0;
        for (int i=0; i<molecules.getAtomCount(); i++) {
            IAtom molecule = molecules.getAtom(i);
            IVector site = getLatticePosition(molecule);
            for (int k = 0; k < site.getD(); k++) {
                work1.setX(k, site.x(k) + newU[j+k]);
            }
            
            atomActionTranslateTo.setDestination(work1);
            atomActionTranslateTo.setAtomPositionDefinition(molecule.getType().getPositionDefinition());
            atomActionTranslateTo.actionPerformed(molecule);
            
            j += coordinateDim/molecules.getAtomCount();

        }

    }

    private static final long serialVersionUID = 1L;
    protected final IVector work1;
    protected final double[] u;
}
