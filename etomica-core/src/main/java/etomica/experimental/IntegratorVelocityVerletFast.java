package etomica.experimental;

import etomica.box.Box;
import etomica.integrator.IntegratorMD;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;

public class IntegratorVelocityVerletFast extends IntegratorMD {
    private final VectorSystem3D positions;
    private final VectorSystem3D forces;

    public IntegratorVelocityVerletFast(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);

        this.positions = new VectorSystem3D(this.box);
        this.forces = new VectorSystem3D(this.box.getLeafList().size());
    }
}
