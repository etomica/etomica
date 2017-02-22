/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;


/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterMoleculeMulti extends MCMoveMolecule {

    protected IVectorRandom[] translationVectors;
    protected int[] constraintMap;
    protected int startMolecule;

    public MCMoveClusterMoleculeMulti(ISimulation sim, ISpace _space) {
    	this(null, sim.getRandom(), _space, 1.0);
    }
    
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in
     * box should be at least one greater than this value (greater
     * because first atom is never moved)
     */
    public MCMoveClusterMoleculeMulti(IPotentialMaster potentialMaster,
            IRandom random, ISpace _space, double stepSize) {
        super(potentialMaster, random, _space, stepSize, Double.POSITIVE_INFINITY);
        setStartMolecule(1);
    }

    public void setBox(IBox p) {
        super.setBox(p);
        translationVectors = new IVectorRandom[box.getMoleculeList().getMoleculeCount()];
        for (int i=0; i<box.getMoleculeList().getMoleculeCount(); i++) {
            translationVectors[i] = (IVectorRandom)space.makeVector();
        }
        if (constraintMap == null) {
            constraintMap = new int[box.getMoleculeList().getMoleculeCount()];
            for (int i=0; i<constraintMap.length; i++) {
                constraintMap[i] = i;
            }
        }
    }
    
    public void setConstraintMap(int[] newConstraintMap) {
        constraintMap = newConstraintMap;
    }

    public void setStartMolecule(int newStartMolecule) {
        startMolecule = newStartMolecule;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
//        if (uOld == 0) {
//            throw new RuntimeException("oops, initial configuration unhappy");
//        }
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=startMolecule; i<moleculeList.getMoleculeCount(); i++) {
            int tv = constraintMap[i];
            if (tv == i) {
                translationVectors[tv].setRandomCube(random);
                translationVectors[tv].TE(stepSize);
            }
            groupTranslationVector.E(translationVectors[tv]);
            moveMoleculeAction.actionPerformed(moleculeList.getMolecule(i));
        }
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
	
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=startMolecule; i<moleculeList.getMoleculeCount(); i++) {
            groupTranslationVector.Ea1Tv1(-1,translationVectors[constraintMap[i]]);
            moveMoleculeAction.actionPerformed(moleculeList.getMolecule(i));
        }
        ((BoxCluster)box).rejectNotify();
        if (((BoxCluster)box).getSampleCluster().value((BoxCluster)box) == 0) {
            throw new RuntimeException("oops oops, reverted to illegal configuration");
        }
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }
    
    public double getB() {
        return 0.0;
    }
    
    public double getA() {
        return uNew/uOld;
    }
	
}
