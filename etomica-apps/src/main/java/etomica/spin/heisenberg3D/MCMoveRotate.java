/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.heisenberg3D;

import etomica.api.*;
import etomica.atom.IAtomOriented;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.IOrientation;
import etomica.space.ISpace;
import etomica.space3d.IOrientation3D;
import etomica.space3d.Orientation3D;

/**
 * Performs a rotation of an atom (not a molecule) that has an orientation coordinate.
 */
public class MCMoveRotate extends MCMoveAtom {
    
    private IOrientation oldOrientation;

    private transient IOrientation iOrientation;

    public MCMoveRotate(IPotentialMaster potentialMaster, IRandom random,
                        ISpace _space) {
        super(potentialMaster, random, _space, Math.PI/2, Math.PI, false);
    }
    
    public void setBox(IBox box) {
        super.setBox(box);
        if (oldOrientation != null) return;
        IAtomList atoms = box.getLeafList();
        if (atoms.getAtomCount() == 0) return;
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        if (atom0.getOrientation() instanceof Orientation3D) {
            oldOrientation = new Orientation3D(space);
        }
        else {
            oldOrientation = space.makeOrientation();
        }
    }

    public boolean doTrial() {
        if(box.getMoleculeList().getMoleculeCount()==0) {return false;}
        atom = atomSource.getAtom();

        energyMeter.setTarget(atom);
        uOld = energyMeter.getDataAsScalar();
        iOrientation = ((IAtomOriented)atom).getOrientation(); 
        oldOrientation.E(iOrientation);  //save old orientation

        //TODO this mcRotate inforce the rotation is around z axis. also the initial
        IVectorMutable dr = space.makeVector();
        dr.setX(2,1);
        ((IOrientation3D)iOrientation).rotateBy((2*random.nextDouble()-1)*stepSize, dr);

        //iOrientation.randomRotation(random, stepSize);


        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public void rejectNotify() {
        iOrientation.E(oldOrientation);
    }
}
