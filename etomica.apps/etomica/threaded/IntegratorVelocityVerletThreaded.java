package etomica.threaded;

import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.util.IRandom;

public class IntegratorVelocityVerletThreaded extends IntegratorVelocityVerlet {

    public IntegratorVelocityVerletThreaded(Simulation sim, int numThreads) {
        this(sim.getPotentialMaster(),sim.getRandom(),
                sim.getDefaults().timeStep,sim.getDefaults().temperature, numThreads);
        // TODO Auto-generated constructor stub
    }

    public IntegratorVelocityVerletThreaded(PotentialMaster potentialMaster,
            IRandom random, double timeStep, double temperature, int numThreads) {
        super(potentialMaster, random, timeStep, temperature);
        
        PotentialCalculationForceSum[] pcfs = new PotentialCalculationForceSum[numThreads];
        for(int i=0; i<numThreads; i++){
            pcfs[i] = new PotentialCalculationForceSum();
        }
        
        forceSum = new PotentialCalculationForceSumThreaded(pcfs);
        // TODO Auto-generated constructor stub
    }

}
