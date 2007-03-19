package etomica.normalmode;

import java.io.Serializable;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * CoordinateDefinition implementation for molecules.  The class takes
 * the first space.D values of u to be real space displacements of the molecule
 * center of mass from its nominal position.  Subclasses should add 
 * additional u values for intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz
 */
public class CoordinateDefinitionMolecule extends CoordinateDefinition implements Serializable {

    public CoordinateDefinitionMolecule(Space space, int orientationDim) {
        super(space.D()+orientationDim);
        this.space = space;
        work1 = space.makeVector();
        atomActionTranslateTo = new AtomActionTranslateTo(space);
    }
    
    public void calcU(Atom molecule, int index, double[] u) {
        IVector pos = molecule.getType().getPositionDefinition().position(molecule);
        for (int i=0; i<pos.getD(); i++) {
            u[i] = pos.x(i) - nominalU[index][i];
        }
    }

    public void initNominalU(Atom molecule, int index) {
        IVector pos = molecule.getType().getPositionDefinition().position(molecule);
        for (int i=0; i<pos.getD(); i++) {
            nominalU[index][i] = pos.x(i);
        }
    }

    public void setToU(Atom molecule, int index, double[] u) {
        for (int i=0; i<space.D(); i++) {
            work1.setX(i, nominalU[index][i] + u[i]);
        }
        atomActionTranslateTo.setDestination(work1);
        atomActionTranslateTo.actionPerformed(molecule);
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final IVector work1;
    protected double[][] nominalU;
    protected final AtomActionTranslateTo atomActionTranslateTo;
}
