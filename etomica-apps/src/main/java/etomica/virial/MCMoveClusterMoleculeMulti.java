/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;


/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterMoleculeMulti extends MCMoveMolecule {

    protected Vector[] translationVectors;
    protected int[] constraintMap;
    protected int startMolecule;

    public MCMoveClusterMoleculeMulti(IRandom random, Space _space) {
        this(null, random, _space, 1.0);
    }

    public MCMoveClusterMoleculeMulti(PotentialMaster potentialMaster,
                                      IRandom random, Space _space, double stepSize) {
        super(potentialMaster, random, _space, stepSize, Double.POSITIVE_INFINITY);
        setStartMolecule(1);
    }

    public void setBox(Box p) {
        super.setBox(p);
        translationVectors = new Vector[box.getMoleculeList().size()];
        for (int i = 0; i<box.getMoleculeList().size(); i++) {
            translationVectors[i] = space.makeVector();
        }
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
        for(int i = startMolecule; i<moleculeList.size(); i++) {
            int tv = constraintMap[i];
            if (tv == i) {
                translationVectors[tv].setRandomCube(random);
                translationVectors[tv].TE(stepSize);
            }
            groupTranslationVector.E(translationVectors[tv]);
            moveMoleculeAction.actionPerformed(moleculeList.get(i));
        }
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
	
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = startMolecule; i<moleculeList.size(); i++) {
            groupTranslationVector.Ea1Tv1(-1,translationVectors[constraintMap[i]]);
            moveMoleculeAction.actionPerformed(moleculeList.get(i));
        }
        ((BoxCluster)box).rejectNotify();
        if (((BoxCluster)box).getSampleCluster().value((BoxCluster)box) == 0) {
            throw new RuntimeException("oops oops, reverted to illegal configuration");
        }
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    public double getChi(double temperature) {
        return uNew/uOld;
    }
	
}
