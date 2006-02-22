package etomica.potential;

import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Potential that wraps a Potential2HardSpherical potential  and fronts 
 * its interface.  Used in cases where it is desired to swap potentials
 * at runtime.  Wrapped potential can be changed while client still uses
 * this instance as its potential. 
 */
public class Potential2HardSphericalWrapper extends Potential2HardSpherical {

    public Potential2HardSphericalWrapper(Space space, Potential2HardSpherical potential) {
        super(space, potential.coordinatePair);
        setPotential(potential);
    }

    public void setPotential(PotentialHard potential) {
        if(potential == null) {
            throw new NullPointerException("Cannot wrap null potential; use P2Ideal");
        }
        if(!(potential instanceof Potential2)) {
            throw new ClassCastException("Can wrap only 2-body potentials");
        }
        if(!(potential instanceof Potential2Spherical)) {
            throw new ClassCastException("Can wrap only spherical potentials");
        }
        wrappedPotential = potential;
    }
    
    public void setPhase(Phase phase) {
        ((Potential2)wrappedPotential).setPhase(phase);
        super.setPhase(phase);
    }
    
    public double getRange() {
        return ((Potential2)wrappedPotential).getRange();
    }
    
    public double collisionTime(AtomSet atoms, double falseTime) {
        return wrappedPotential.collisionTime(atoms,falseTime);
    }
    
    public double lastCollisionVirial() {
        return wrappedPotential.lastCollisionVirial();
    }
    
    public Tensor lastCollisionVirialTensor() {
        return wrappedPotential.lastCollisionVirialTensor();
    }
    
    public void bump(AtomSet atoms, double falseTime) {
        wrappedPotential.bump(atoms,falseTime);
    }
    
    public double energyChange() {
        return wrappedPotential.energyChange();
    }

    public double u(double r2) {
        return ((Potential2Spherical)wrappedPotential).u(r2);
    }
    
    private PotentialHard wrappedPotential;

}
