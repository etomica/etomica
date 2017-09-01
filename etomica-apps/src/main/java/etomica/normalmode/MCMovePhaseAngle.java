/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.iterator.AtomIterator;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMovePhaseAngle extends MCMoveBoxStep {
    
    public MCMovePhaseAngle(Space space, NormalModesVariable normalModes, CoordinateDefinition coordinateDefinition, IRandom random) {
        super(null);
        this.space = space;
        this.normalModes = normalModes;
        this.coordinateDefinition = coordinateDefinition;
        this.random = random;
        meterJacobian = new MeterJacobian(coordinateDefinition, normalModes);
        init();
        setStepSize(1);
        setStepSizeMax(Math.PI/2);
    }

    protected void init() {
        oldPhaseAngles = new double[normalModes.getWaveVectors().length];
    }
    
    public boolean doTrial() {
        Vector[] waveVectors = normalModes.getWaveVectors();
        double[] phaseAngles = normalModes.getPhaseAngles();
        if (phaseAngles.length != oldPhaseAngles.length) {
            init();
        }
        for (int i=0; i<waveVectors.length; i++) {
            oldPhaseAngles[i] = phaseAngles[i];
        }
        
        jOld = meterJacobian.getDataAsScalar();
        
        int mode = random.nextInt(waveVectors.length-1)+1;
        phaseAngles[mode] += 2.0*stepSize*(random.nextDouble()-0.5);
        if (phaseAngles[mode] < 0) {
            phaseAngles[mode] += 2.0*Math.PI;
        }
        else if (phaseAngles[mode] > 2.0*Math.PI) {
            phaseAngles[mode] -= 2.0*Math.PI;
        }
            
//            System.out.println("changing mode # "+mode);
//            System.out.println("theta "+oldPhaseAngles[mode]+" => "+phaseAngles[mode]);
        
        jNew = meterJacobian.getDataAsScalar();
        
        System.out.println("jOld "+jOld+"    jNew "+jNew);
        
        return true;
    }

    public double getChi(double temperature) {
        return jNew/jOld;
    }

    public void acceptNotify() {
//        System.out.println("accepted");
    }

    public void rejectNotify() {
//        System.out.println("rejected");
        double[] phaseAngles = normalModes.getPhaseAngles();
        for (int i=0; i<phaseAngles.length; i++) {
            phaseAngles[i] = oldPhaseAngles[i];
        }
    }

    public AtomIterator affectedAtoms() {
        return null;
    }

    public double energyChange() {
        return 0;
    }

    protected final Space space;
    protected final NormalModesVariable normalModes;
    protected double[] oldPhaseAngles;
    protected final CoordinateDefinition coordinateDefinition;
    protected final IRandom random;
    protected MeterJacobian meterJacobian;
    protected double jOld, jNew;
}
