/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveWV extends MCMoveBoxStep {
    
    public MCMoveWV(Space space, NormalModesVariable normalModes, CoordinateDefinition coordinateDefinition, IRandom random) {
        super(null);
        this.space = space;
        this.normalModes = normalModes;
        this.coordinateDefinition = coordinateDefinition;
        this.random = random;
        meterJacobian = new MeterJacobian(coordinateDefinition, normalModes);
        init();
    }

    protected void init() {
        oldEigenVectors = new double[normalModes.getEigenVectors().length][coordinateDefinition.getCoordinateDim()];
        oldWaveVectors = new Vector[normalModes.getWaveVectors().length];
        for (int i=0; i<oldWaveVectors.length; i++) {
            oldWaveVectors[i] = space.makeVector();
        }
    }

    public void setBox(Box newBox) {
        super.setBox(newBox);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        double L = box.getBoundary().getBoxSize().getX(0);
        double minWV = 2. * Math.PI / L;
        double maxWV = 2. * Math.PI * (nMolecules/2) / L;
        setStepSize(0.2*(maxWV-minWV));
        setStepSizeMax(0.25*(maxWV-minWV));
    }

    public boolean doTrial() {
        Vector[] waveVectors = normalModes.getWaveVectors();
        double[][] eigenVectors = normalModes.getEigenVectors();
        if (waveVectors.length != oldWaveVectors.length) {
            init();
        }
        for (int i=0; i<waveVectors.length; i++) {
            oldWaveVectors[i].E(waveVectors[i]);
            for (int j=0; j<coordinateDefinition.getCoordinateDim(); j++) {
                oldEigenVectors[i][j] = eigenVectors[i][j];
            }
        }
        
        jOld = meterJacobian.getDataAsScalar();
        if (jOld < 0.001) {
            System.out.println("wv:");
            for (int i=0; i<waveVectors.length; i++) {
                System.out.println(waveVectors[i]+ " "+ normalModes.getPhaseAngles()[i]);
            }
        }
        
        int mode = random.nextInt(waveVectors.length-1)+1;
        if (waveVectors[0].getD() == 1) {
            // density is 1
            double wv = waveVectors[mode].getX(0);
            wv += stepSize*2.0*(random.nextDouble()-0.5);
            int nMolecules = box.getMoleculeList().getMoleculeCount();
            double L = box.getBoundary().getBoxSize().getX(0);
            double minWV = 2. * Math.PI / L;
            double maxWV = 2. * Math.PI * (nMolecules/2) / L;
            if (wv < minWV) wv += (maxWV-minWV);
            else if (wv > maxWV) wv -= (maxWV-minWV);
            waveVectors[mode].E(wv);
            // skip eigenvector
//            phaseAngles[mode] = 2.0*Math.PI*random.nextDouble();
//            System.out.println("changing mode # "+mode);
//            System.out.println("WV "+oldWaveVectors[mode].getX(0)+" => "+waveVectors[mode].getX(0));
//            System.out.println("theta "+oldPhaseAngles[mode]+" => "+phaseAngles[mode]);
        }
        
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
        Vector[] waveVectors = normalModes.getWaveVectors();
        double[][] eigenVectors = normalModes.getEigenVectors();
        for (int i=0; i<waveVectors.length; i++) {
            waveVectors[i].E(oldWaveVectors[i]);
            for (int j=0; j<coordinateDefinition.getCoordinateDim(); j++) {
                eigenVectors[i][j] = oldEigenVectors[i][j];
            }
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
    protected Vector[] oldWaveVectors;
    protected double[][] oldEigenVectors;
    protected final CoordinateDefinition coordinateDefinition;
    protected final IRandom random;
    protected MeterJacobian meterJacobian;
    protected double jOld, jNew;
}
