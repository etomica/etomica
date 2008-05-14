package etomica.potential;

import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Dimension;

/**
 * Hard potential class that wraps another hard potential.
 *
 * @author Andrew Schultz
 */
 public class P2HardWrapper implements PotentialHard {

    public P2HardWrapper(ISpace space, PotentialHard potential) {
        this.space = space;
        wrappedPotential = potential;
    }

    public double energy(IAtomSet atoms) {
        return wrappedPotential.energy(atoms);
    }

    public int nBody() {
        return wrappedPotential.nBody();
    }

    public void setBox(IBox box) {
        wrappedPotential.setBox(box);
    }
    
    public PotentialHard getWrappedPotential() {
        return wrappedPotential;
    }

    public void setWrappedPotential(PotentialHard newWrappedPotential) {
        wrappedPotential = newWrappedPotential;
    }

    public double getRange() {
        return wrappedPotential.getRange();
    }
    
    public void bump(IAtomSet atoms, double falseTime) {
        wrappedPotential.bump(atoms, falseTime);
    }

    public double collisionTime(IAtomSet atoms, double falseTime) {
        return wrappedPotential.collisionTime(atoms, falseTime);
    }

    public double energyChange() {
        return wrappedPotential.energyChange();
    }

    public double lastCollisionVirial() {
        return wrappedPotential.lastCollisionVirial();
    }

    public Tensor lastCollisionVirialTensor() {
        return wrappedPotential.lastCollisionVirialTensor();
    }

    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected PotentialHard wrappedPotential;
}
