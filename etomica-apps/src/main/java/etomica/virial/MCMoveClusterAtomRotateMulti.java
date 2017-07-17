/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.util.random.IRandom;
import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.Orientation3D;

/**
 * Extension of MCMoveAtom that does trial in which several atom orientations are
 * perturbed.  However, orientation of first atom is never altered.  
 */
public class MCMoveClusterAtomRotateMulti extends MCMoveAtom {

    public MCMoveClusterAtomRotateMulti(IRandom random, Space _space) {
        super(random, null, _space, 1.0, Math.PI, false);
        setStepSize(1.2);
	}

    public void setBox(Box box) {
        super.setBox(box);
        if (oldOrientations != null) return;
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.getAtomCount();
        oldOrientations = new IOrientation[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            if (((IAtomOriented)atoms.getAtom(i)).getOrientation() instanceof Orientation3D) {
                oldOrientations[i] = new Orientation3D(space);
            }
            else {
                oldOrientations[i] = space.makeOrientation();
            }
        }
    }
    
	//note that total energy is calculated
	public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.getAtomCount();
        for(int i=0; i<nAtoms; i++) {
            IAtomOriented a = (IAtomOriented)atoms.getAtom(i);
            oldOrientations[i].E(a.getOrientation());
            a.getOrientation().randomRotation(random, stepSize);
        }
        if (random.nextInt(100) == 0) {
            for(int i=0; i<nAtoms; i++) {
                IAtomOriented a = (IAtomOriented)atoms.getAtom(i);
                a.getOrientation().getDirection().normalize();
            }
        }
		((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
		return true;
	}
	
    public double getA() {
        return uNew/uOld;
    }

    public double getB() {
    	return 0.0;
    }

    public void rejectNotify() {
        IAtomList atoms = box.getLeafList();
        int nAtoms = atoms.getAtomCount();
        for(int i=0; i<nAtoms; i++) {
            IAtomOriented a = (IAtomOriented)atoms.getAtom(i);
            a.getOrientation().E(oldOrientations[i]);
        }
    	((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    private IOrientation[] oldOrientations;
}
