package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.api.IVector;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * 2-body soft Potential class for use between two monatomic molecules.
 *
 * @author Andrew Schultz
 */
 public class P2SoftMoleculeMonatomic extends P2MoleculeMonatomic implements
        Potential2Soft {

    public P2SoftMoleculeMonatomic(Space space, IPotential potential) {
        super(space, potential);
    }

    public double hyperVirial(IAtomSet pair) {
        leafAtoms.atom0 = ((IMolecule)pair.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)pair.getAtom(1)).getChildList().getAtom(0);
        return wrappedPotential.energy(leafAtoms);
    }

    public double integral(double rc) {
        return ((Potential2Soft)wrappedPotential).integral(rc);
    }

    public double u(double r2) {
        return ((Potential2Soft)wrappedPotential).u(r2);
    }

    public IVector[] gradient(IAtomSet atoms) {
        leafAtoms.atom0 = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)atoms.getAtom(1)).getChildList().getAtom(0);
        return ((PotentialSoft)wrappedPotential).gradient(leafAtoms);
    }

    public IVector[] gradient(IAtomSet atoms, Tensor pressureTensor) {
        leafAtoms.atom0 = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)atoms.getAtom(1)).getChildList().getAtom(0);
        return ((PotentialSoft)wrappedPotential).gradient(leafAtoms, pressureTensor);
    }

    public double virial(IAtomSet atoms) {
        leafAtoms.atom0 = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)atoms.getAtom(1)).getChildList().getAtom(0);
        return ((PotentialSoft)wrappedPotential).virial(leafAtoms);
    }

    private static final long serialVersionUID = 1L;
}
