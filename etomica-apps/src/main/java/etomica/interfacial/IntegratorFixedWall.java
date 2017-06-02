package etomica.interfacial;

import etomica.potential.PotentialMaster;
import etomica.api.IRandom;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.space.Space;

public class IntegratorFixedWall extends IntegratorVelocityVerlet {

    protected FixedWall fixedWall;
    
    public IntegratorFixedWall(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Space space) {
        super(potentialMaster, random, timeStep, temperature, space);
    }
    
    public void setFixedWall(FixedWall fixedWall) {
        this.fixedWall = fixedWall;
        getEventManager().addListener(fixedWall);
    }
    
    public void randomizeMomenta() {
        super.randomizeMomenta();
        if (fixedWall != null) fixedWall.integratorInitialized(null);
    }
}
