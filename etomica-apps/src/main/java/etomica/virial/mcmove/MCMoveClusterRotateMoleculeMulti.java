/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.action.MoleculeAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.RotationTensor;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * MCMove for use in a Mayer sampling simulation that rotates all molecules in
 * a Box except the first molecule, which is never moved.  The angle of
 * rotation is the step size and can be tuned for some acceptance rate.
 */
public class MCMoveClusterRotateMoleculeMulti extends MCMoveBoxStep {

    protected double wOld, wNew;
    protected final IRandom random;
    protected int[] constraintMap;
    protected int[] rotationAxis;
    protected double[] theta;
    private int trialCount, relaxInterval = 100;
    private MoleculeAction relaxAction;
    protected final Vector r0;
    protected final RotationTensor rotationTensor;
    protected boolean doLattice;
    protected int[] rotationCenters;

    public MCMoveClusterRotateMoleculeMulti(IRandom random, Box box) {
        super();
        this.random = random;
        setBox(box);
        r0 = box.getSpace().makeVector();
        rotationTensor = box.getSpace().makeRotationTensor();
        setStepSizeMax(Math.PI);
    }

    public void setDoLattice(boolean doLattice) {
        this.doLattice = doLattice;
    }

    public void setBox(Box p) {
        super.setBox(p);
        IMoleculeList moleculeList = box.getMoleculeList();
        rotationAxis = new int[moleculeList.size()];
        theta = new double[rotationAxis.length];
        rotationCenters = new int[rotationAxis.length];
        if (constraintMap == null) {
            constraintMap = new int[box.getMoleculeList().size()];
            for (int i=0; i<constraintMap.length; i++) {
                constraintMap[i] = i;
            }
        }
    }

    @Override
    public double energyChange() {
        return 0;
    }

    public void setConstraintMap(int[] newConstraintMap) {
        constraintMap = newConstraintMap;
    }

    protected void doTransform(IMolecule molecule) {
        IAtomList childList = molecule.getChildList();
        for (IAtom a : childList) {
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
            r.PE(box.getBoundary().centralImage(r));
        }
    }

    public boolean doTrial() {
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
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
            IMolecule molecule = moleculeList.get(i);
            if (doLattice) {
                if (i==0) continue;
                rotationCenters[i] = random.nextInt(box.getLeafList().size());
                r0.E(box.getLeafList().get(rotationCenters[i]).getPosition());
            }
            else {
                r0.E(CenterOfMass.position(box, molecule));
            }

            int j = constraintMap[i];
            if (j == i) {
                if (doLattice) {
                    theta[j] = Math.PI/2 * (random.nextInt(2)*2-1);
                }
                else {
                    theta[j] = (2 * random.nextDouble() - 1.0) * stepSize;
                }
                rotationAxis[j] = random.nextInt(3);
            }
            rotationTensor.setAxial(rotationAxis[j],theta[j]);

            doTransform(molecule);

            if (doRelax && relaxAction != null) {
                relaxAction.actionPerformed(molecule);
            }

            if (doLattice && (i==0 || (constraintMap != null && constraintMap[i]==0))) {
                Vector shift = box.getSpace().makeVector();
                shift.E(CenterOfMass.position(box, molecule));
                shift.TE(-1);
                for (int k=0; k<shift.getD(); k++) {
                    shift.setX(k, Math.floor(shift.getX(k)));
                }
                for (IAtom aa : molecule.getChildList()) {
                    aa.getPosition().PE(shift);
                }
            }

        }

        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    public double getChi(double temperature) {
        return wNew/wOld;
    }
    
    public void acceptNotify() {
//        if (uNew == 0) {
//            throw new RuntimeException("oops, accepted illegal configuration");
//        }
        ((BoxCluster)box).acceptNotify();
    }
    
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for (int i = 0; i<moleculeList.size(); i++) {
            IMolecule molecule = moleculeList.get(i);
            if (doLattice) {
                if (i==0) continue;
                r0.E(box.getLeafList().get(rotationCenters[i]).getPosition());
            }
            else {
                r0.E(CenterOfMass.position(box, molecule));
            }
            int j = constraintMap[i];
            rotationTensor.setAxial(rotationAxis[j],-theta[j]);

            doTransform(molecule);

            if (doLattice && (i==0 || (constraintMap != null && constraintMap[i]==0))) {
                Vector shift = box.getSpace().makeVector();
                shift.E(CenterOfMass.position(box, molecule));
                shift.TE(-1);
                for (int k=0; k<shift.getD(); k++) {
                    shift.setX(k, Math.floor(shift.getX(k)));
                }
                for (IAtom aa : molecule.getChildList()) {
                    aa.getPosition().PE(shift);
                }
            }
        }
        ((BoxCluster)box).rejectNotify();
        if (((BoxCluster)box).getSampleCluster().value((BoxCluster)box) == 0) {
            throw new RuntimeException("oops oops, reverted to illegal configuration");
        }
    }
    
    public void setRelaxAction(MoleculeAction action) {
        relaxAction = action;
    }
}
