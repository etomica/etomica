package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.atom.AtomPair;
import etomica.atom.AtomSetSinglet;
import etomica.space.Space;

/**
 * 2-body Potential class for use between two monatomic molecules.  The
 * potential can be used instead of a PotentialGroup (which has very high
 * overhead).
 *
 * @author Andrew Schultz
 */
public class P1MoleculeMonatomic extends Potential {

    public P1MoleculeMonatomic(Space space, IPotential potential) {
        super(1, space);
        wrappedPotential = potential;
        leafAtomSet = new AtomSetSinglet();
    }
    
    public double energy(IAtomSet atoms) {
        leafAtomSet.atom = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        return wrappedPotential.energy(leafAtomSet);
    }

    public double getRange() {
        return wrappedPotential.getRange();
    }

    public void setBox(IBox box) {
        wrappedPotential.setBox(box);
    }
    
    public IPotential getWrappedPotential() {
        return wrappedPotential;
    }

    public void setWrappedPotential(IPotential newWrappedPotential) {
        wrappedPotential = newWrappedPotential;
    }

    private static final long serialVersionUID = 1L;
    protected final AtomSetSinglet leafAtomSet;
    protected IPotential wrappedPotential;
}
