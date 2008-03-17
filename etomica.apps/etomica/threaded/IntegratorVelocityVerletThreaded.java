package etomica.threaded;

import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.space.Space;


public class IntegratorVelocityVerletThreaded extends IntegratorVelocityVerlet {

    private static final long serialVersionUID = 1L;

    public IntegratorVelocityVerletThreaded(IRandom _random, Space _space, IPotentialMaster potentialMaster, int numThreads) {
        this(potentialMaster,_random, 0.05, 1.0, numThreads, _space);
    }

    public IntegratorVelocityVerletThreaded(IPotentialMaster potentialMaster,
            IRandom random, double timeStep, double temperature, int numThreads, Space _space) {
        super(potentialMaster, random, timeStep, temperature, _space);
        
        PotentialCalculationForceSum[] pcfs = new PotentialCalculationForceSum[numThreads];
        for(int i=0; i<numThreads; i++){
            pcfs[i] = new PotentialCalculationForceSum();
        }
        
        forceSum = new PotentialCalculationForceSumThreaded(pcfs, _space);
    }

}
