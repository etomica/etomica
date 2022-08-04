/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.IOrientation;
import etomica.space3d.Orientation3D;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

/**
 * Extension of MCMoveAtom that does trial in which several atom orientations are
 * perturbed.  However, orientation of first atom is never altered.  
 */
public class MCMoveClusterAtomRotateMulti extends MCMoveBoxStep {

    protected double wOld, wNew;
    protected final IRandom random;

    public MCMoveClusterAtomRotateMulti(IRandom random, Box box) {
        super();
        this.random = random;
        setStepSize(1.2);
        setStepSizeMax(Math.PI);
        setBox(box);
	}

    public void setBox(Box box) {
        super.setBox(box);
        if (oldOrientations != null) return;
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.size();
        oldOrientations = new IOrientation[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            if (((IAtomOriented)atoms.get(i)).getOrientation() instanceof Orientation3D) {
                oldOrientations[i] = new Orientation3D(box.getSpace());
            }
            else {
                oldOrientations[i] = box.getSpace().makeOrientation();
            }
        }
    }

    @Override
    public double energyChange() {
        return 0;
    }

    //note that total energy is calculated
	public boolean doTrial() {
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.size();
        for(int i=0; i<nAtoms; i++) {
            IAtomOriented a = (IAtomOriented)atoms.get(i);
            oldOrientations[i].E(a.getOrientation());
            a.getOrientation().randomRotation(random, stepSize);
        }
        if (random.nextInt(100) == 0) {
            for(int i=0; i<nAtoms; i++) {
                IAtomOriented a = (IAtomOriented)atoms.get(i);
                a.getOrientation().getDirection().normalize();
            }
        }
		((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}

    public double getChi(double temperature) {
        return wNew/wOld;
    }

    public void rejectNotify() {
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.size();
        for(int i=0; i<nAtoms; i++) {
            IAtomOriented a = (IAtomOriented)atoms.get(i);
            a.getOrientation().E(oldOrientations[i]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    private IOrientation[] oldOrientations;
}
