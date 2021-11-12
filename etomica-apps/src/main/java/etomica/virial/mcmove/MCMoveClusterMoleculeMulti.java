/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;


/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterMoleculeMulti extends MCMoveBoxStep {

    protected double uOld, uNew;
    protected final IRandom random;
    protected Vector[] translationVectors;
    protected int[] constraintMap;
    protected int startMolecule;

    public MCMoveClusterMoleculeMulti(IRandom random, Box box) {
        super();
        this.random = random;
        setBox(box);
        setStartMolecule(1);
    }

    public void setBox(Box p) {
        super.setBox(p);
        translationVectors = new Vector[box.getMoleculeList().size()];
        for (int i = 0; i<box.getMoleculeList().size(); i++) {
            translationVectors[i] = box.getSpace().makeVector();
        }
        if (constraintMap == null) {
            constraintMap = new int[box.getMoleculeList().size()];
            for (int i=0; i<constraintMap.length; i++) {
                constraintMap[i] = i;
            }
        }
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
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
            moleculeList.get(i).getChildList().forEach(atom -> {
                atom.getPosition().PE(translationVectors[tv]);
            });
        }
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
	
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = startMolecule; i<moleculeList.size(); i++) {
            Vector v = translationVectors[constraintMap[i]];
            moleculeList.get(i).getChildList().forEach(atom -> {
                atom.getPosition().ME(v);
            });
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
