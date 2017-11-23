/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * @author kofke
 *
 * Extension of MCMoveAtom that does trial in which several atom positions are
 * perturbed.  However, position of first atom is never altered.  
 */
public class MCMoveClusterAtomMulti extends MCMoveAtom {

    public MCMoveClusterAtomMulti(IRandom random, Space _space) {
        super(random, null, _space);
        setStepSize(1.2);
	}
	
    public void setBox(Box p) {
        super.setBox(p);
        if (translationVectors == null) {
            translationVectors = new Vector[box.getLeafList().getAtomCount()-1];
            for (int i=0; i<translationVectors.length; i++) {
                translationVectors[i] = space.makeVector();
            }
        }
    }
    
    public void setStartAtom(int newStartAtom) {
        startAtom = newStartAtom;
    }
    
    public int getStartAtom() {
        return startAtom;
    }
    
	public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IAtomList leafAtoms = box.getLeafList();
        for(int i=startAtom; i<leafAtoms.getAtomCount(); i++) {
            translationVectors[i-1].setRandomCube(random);
            translationVectors[i-1].TE(stepSize);
            leafAtoms.getAtom(i).getPosition().PE(translationVectors[i-1]);
        }
		((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}

    public double getChi(double temperature) {
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public void rejectNotify() {
        IAtomList leafAtoms = box.getLeafList();
        for(int i=startAtom; i<leafAtoms.getAtomCount(); i++) {
            leafAtoms.getAtom(i).getPosition().ME(translationVectors[i-1]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    protected Vector[] translationVectors;
    protected int startAtom = 1;
}
