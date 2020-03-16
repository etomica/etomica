/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.MoleculeAction;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * MCMove for use in a Mayer sampling simulation that rotates all molecules in
 * a Box except the first molecule, which is never moved.  The angle of
 * rotation is the step size and can be tuned for some acceptance rate.
 */
public class MCMoveClusterRotateMoleculeMulti extends MCMoveRotateMolecule3D {

    public MCMoveClusterRotateMoleculeMulti(IRandom random, Space _space) {
        super(null, random, _space);
    }

    public void setBox(Box p) {
        super.setBox(p);
        IMoleculeList moleculeList = box.getMoleculeList();
        rotationAxis = new int[moleculeList.size()];
        theta = new double[rotationAxis.length];
        if (constraintMap == null) {
            constraintMap = new int[box.getMoleculeList().size()];
            for (int i=0; i<constraintMap.length; i++) {
                constraintMap[i] = i;
            }
        }
    }

    public void setConstraintMap(int[] newConstraintMap) {
        constraintMap = newConstraintMap;
    }

    public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
//        if (uOld == 0) {
//            throw new RuntimeException("oops, illegal initial configuration");
//        }
        boolean doRelax = false;
        if (trialCount-- == 0) {
            doRelax = true;
            trialCount = relaxInterval;
        }
        IMoleculeList moleculeList = box.getMoleculeList();
        for (int i = 0; i<moleculeList.size(); i++) {
            molecule = moleculeList.get(i);
            r0.E(positionDefinition.position(molecule));

            int j = constraintMap[i];
            if (j == i) {
                theta[j] = (2*random.nextDouble() - 1.0)*stepSize;
                rotationAxis[j] = random.nextInt(3);
            }

            rotationTensor.setAxial(rotationAxis[j],theta[j]);

            doTransform();

            if (doRelax && relaxAction != null) {
                relaxAction.actionPerformed(molecule);
            }
        }

        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    public double getChi(double temperature) {
        return uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
//        if (uNew == 0) {
//            throw new RuntimeException("oops, accepted illegal configuration");
//        }
        ((BoxCluster)box).acceptNotify();
    }
    
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for (int i = 0; i<moleculeList.size(); i++) {
            molecule = moleculeList.get(i);
            r0.E(positionDefinition.position(molecule));
            int j = constraintMap[i];
            rotationTensor.setAxial(rotationAxis[j],-theta[j]);

            doTransform();
        }
        ((BoxCluster)box).rejectNotify();
        if (((BoxCluster)box).getSampleCluster().value((BoxCluster)box) == 0) {
            throw new RuntimeException("oops oops, reverted to illegal configuration");
        }
    }
    
    public void setRelaxAction(MoleculeAction action) {
        relaxAction = action;
    }
    
    private static final long serialVersionUID = 1L;
    protected int[] constraintMap;
    protected int[] rotationAxis;
    protected double[] theta;
    private int trialCount, relaxInterval = 100;
    private MoleculeAction relaxAction;
}
