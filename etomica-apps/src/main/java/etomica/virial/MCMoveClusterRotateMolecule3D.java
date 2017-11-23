/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.action.MoleculeAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.util.random.IRandom;

/**
 * MC move for Mayer Sampling that rotates a single molecule
 */
public class MCMoveClusterRotateMolecule3D extends MCMoveRotateMolecule3D {

    public MCMoveClusterRotateMolecule3D(PotentialMaster potentialMaster,
            IRandom random, Space _space) {
        super(potentialMaster, random, _space);
    }

    public boolean doTrial() {
        molecule = moleculeSource.getMolecule();
        while (molecule.getIndex() == 0) {
            molecule = moleculeSource.getMolecule();
        }
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        
        double dTheta = (2*random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(random.nextInt(3),dTheta);

        IAtomList leafAtoms = molecule.getChildList();
        IAtom first = leafAtoms.getAtom(0);
        r0.E(first.getPosition());
        doTransform();

        if (trialCount-- == 0) {
            relaxAction.actionPerformed(molecule);
            trialCount = relaxInterval;
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
        ((BoxCluster)box).acceptNotify();
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        ((BoxCluster)box).rejectNotify();
    }
    
    public void setRelaxAction(MoleculeAction action) {
        relaxAction = action;
    }
    
    private static final long serialVersionUID = 1L;
    protected int trialCount, relaxInterval = 100;
    protected MoleculeAction relaxAction;
}
