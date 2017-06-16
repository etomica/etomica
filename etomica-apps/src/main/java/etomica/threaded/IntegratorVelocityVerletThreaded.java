/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.threaded;

import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.space.Space;


public class IntegratorVelocityVerletThreaded extends IntegratorVelocityVerlet {

    private static final long serialVersionUID = 1L;

    public IntegratorVelocityVerletThreaded(IRandom _random, Space _space, PotentialMaster potentialMaster, int numThreads) {
        this(potentialMaster,_random, 0.05, 1.0, numThreads, _space);
    }

    public IntegratorVelocityVerletThreaded(PotentialMaster potentialMaster,
            IRandom random, double timeStep, double temperature, int numThreads, Space _space) {
        super(potentialMaster, random, timeStep, temperature, _space);
        
        PotentialCalculationForceSum[] pcfs = new PotentialCalculationForceSum[numThreads];
        for(int i=0; i<numThreads; i++){
            pcfs[i] = new PotentialCalculationForceSum();
        }
        
        forceSum = new PotentialCalculationForceSumThreaded(pcfs, _space);
    }

}
