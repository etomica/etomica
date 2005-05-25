package etomica.potential;

import etomica.AtomSet;
import etomica.Space;
import etomica.space.Tensor;

/**
 * Potential that wraps a Potential2HardSpherical potential  and fronts 
 * its interface.  Used in cases where it is desired to swap potentials
 * at runtime.  Wrapped potential can be changed while client still uses
 * this instance as its potential. 
 */
public class Potential2HardSphericalWrapper extends Potential2HardSpherical {

    public Potential2HardSphericalWrapper(Space space, Potential2HardSpherical potential) {
        super(space);
        setPotential(potential);
    }

    public void setPotential(Potential2HardSpherical potential) {
        wrappedPotential = potential;
    }
    
    public double getRange() {
        return wrappedPotential.getRange();
    }
    
    public double collisionTime(AtomSet atoms, double falseTime) {
        if (wrappedPotential==null) return Double.POSITIVE_INFINITY;
        return wrappedPotential.collisionTime(atoms,falseTime);
    }
    
    public double lastCollisionVirial() {
        if (wrappedPotential==null) throw new RuntimeException("can't have virial with null potential");;
        return wrappedPotential.lastCollisionVirial();
    }
    
    public Tensor lastCollisionVirialTensor() {
        if (wrappedPotential==null) throw new RuntimeException("can't have virial with null potential");
        return wrappedPotential.lastCollisionVirialTensor();
    }
    
    public void bump(AtomSet atoms, double falseTime) {
        if (wrappedPotential==null) throw new RuntimeException("can't bump with null potential");
        wrappedPotential.bump(atoms,falseTime);
    }
    
    public double energyChange() {
        if (wrappedPotential==null) throw new RuntimeException("can't have energy change with null potential");
        return wrappedPotential.energyChange();
    }

    public double u(double r2) {
        if (wrappedPotential==null) return 0;
        return wrappedPotential.u(r2);
    }
    
    private Potential2HardSpherical wrappedPotential;

    
}
