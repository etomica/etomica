/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.util.random.IRandom;
import etomica.atom.IAtomOriented;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.Orientation3D;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */
public class MCMoveRotate extends MCMoveAtom {
    
    private IOrientation oldOrientation;

    private transient IOrientation iOrientation;

    public MCMoveRotate(PotentialMaster potentialMaster, IRandom random,
                        Space _space) {
        super(random, potentialMaster, _space, Math.PI/2, Math.PI, false);
    }
    
    public void setBox(Box box) {
        super.setBox(box);
        if (oldOrientation != null) return;
        IAtomList atoms = box.getLeafList();
        if (atoms.size() == 0) return;
        IAtomOriented atom0 = (IAtomOriented)atoms.get(0);
        if (atom0.getOrientation() instanceof Orientation3D) {
            oldOrientation = new Orientation3D(space);
        }
        else {
            oldOrientation = space.makeOrientation();
        }
    }

    public boolean doTrial() {
        if(box.getMoleculeList().size()==0) {return false;}
        atom = atomSource.getAtom();

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        iOrientation = ((IAtomOriented)atom).getOrientation(); 
        oldOrientation.E(iOrientation);  //save old orientation
        iOrientation.randomRotation(random, stepSize);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public void rejectNotify() {
        iOrientation.E(oldOrientation);
    }
}
