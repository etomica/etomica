package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IMolecule;
import etomica.api.IPotential;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * 2-body hard Potential class for use between two monatomic molecules.
 *
 * @author Andrew Schultz
 */
 public class P2HardMoleculeMonatomic extends P2MoleculeMonatomic implements
        PotentialHard {

    public P2HardMoleculeMonatomic(Space space, IPotential potential) {
        super(space, potential);
    }

    public void bump(IAtomSet atoms, double falseTime) {
        leafAtoms.atom0 = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)atoms.getAtom(1)).getChildList().getAtom(0);
        ((PotentialHard)wrappedPotential).bump(leafAtoms, falseTime);
    }

    public double collisionTime(IAtomSet atoms, double falseTime) {
        leafAtoms.atom0 = ((IMolecule)atoms.getAtom(0)).getChildList().getAtom(0);
        leafAtoms.atom1 = ((IMolecule)atoms.getAtom(1)).getChildList().getAtom(0);
        return ((PotentialHard)wrappedPotential).collisionTime(leafAtoms, falseTime);
    }

    public double energyChange() {
        return ((PotentialHard)wrappedPotential).energyChange();
    }

    public double lastCollisionVirial() {
        return ((PotentialHard)wrappedPotential).lastCollisionVirial();
    }

    public Tensor lastCollisionVirialTensor() {
        return ((PotentialHard)wrappedPotential).lastCollisionVirialTensor();
    }

    private static final long serialVersionUID = 1L;
}
