/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
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
public class MCMoveClusterMoleculeFlipSide extends MCMoveBoxStep {

    protected double uOld, uNew;
    protected final IRandom random;
    protected Vector translationVector;
    protected boolean doLattice;
    protected IMolecule molecule;

    public MCMoveClusterMoleculeFlipSide(IRandom random, Box box) {
        super();
        this.random = random;
        setBox(box);
    }

    public void setBox(Box p) {
        super.setBox(p);
        translationVector = box.getSpace().makeVector();
    }


    @Override
    public double energyChange() {
        return 0;
    }



    //note that total energy is calculated
    public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
//        if (uOld == 0) {
//            throw new RuntimeException("oops, initial configuration unhappy");
//        }
        IMoleculeList moleculeList = box.getMoleculeList();
        if(moleculeList.size() == 4){
            molecule = moleculeList.get(random.nextInt(2) +1);

        }
        else{
            molecule = moleculeList.get(1);
        }
        translationVector.Ea1Tv1(-2, CenterOfMass.position(box, molecule));
//        translationVector.setX(1, 0);
//        translationVector.setX(2, 0);
        molecule.getChildList().forEach(atom -> {
            atom.getPosition().PE(translationVector);
        });


        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
	
    public void rejectNotify() {

        molecule.getChildList().forEach(atom -> {
            atom.getPosition().ME(translationVector);
        });

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
