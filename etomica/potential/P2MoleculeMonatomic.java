package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.atom.AtomPair;
import etomica.space.Space;

/**
 * 2-body Potential class for use between two monatomic molecules.  The
 * potential can be used instead of a PotentialGroup (which has very high
 * overhead).
 *
 * @author Andrew Schultz
 */
public class P2MoleculeMonatomic extends Potential {

    public P2MoleculeMonatomic(Space space, IPotential potential) {
        super(2, space);
        wrappedPotential = potential;
        leafAtoms = new AtomPair();
    }
    
    public double energy(IAtomSet atoms) {
        leafAtoms.atom0 = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)atoms.getAtom(1)).getChildList().getAtom(0);
        return wrappedPotential.energy(leafAtoms);
    }

    public double getRange() {
        return wrappedPotential.getRange();
    }

    public void setBox(IBox box) {
        wrappedPotential.setBox(box);
    }

    private static final long serialVersionUID = 1L;
    protected final AtomPair leafAtoms;
    protected final IPotential wrappedPotential;
}
