package etomica.potential;

import etomica.AtomPair;
import etomica.Space;
import etomica.space.Tensor;

public class Potential2HardSphericalWrapper extends Potential2HardSpherical {

    public Potential2HardSphericalWrapper(Space space, Potential2HardSpherical potential) {
        super(space);
        setPotential(potential);
    }

    public void setPotential(Potential2HardSpherical potential) {
        wrappedPotential = potential;
    }
    
    public double collisionTime(AtomPair atoms, double falseTime) {
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
    
    public void bump(AtomPair atoms, double falseTime) {
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
