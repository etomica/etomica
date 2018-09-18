/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * An MC Move for cluster simulations that performs torsion moves on acetic acid.
 * The move is performed on all molecules in the Box.  
 * The move needs the torsion potential in order to choose appropriate
 * torsion angles.
 * Only the hydrogen site is rotating by 180 degree. 
 * 
 * @author Hye Min Kim
 */
public class MCMoveTorsionAceticAcid extends MCMoveMolecule {
   
    public MCMoveTorsionAceticAcid(PotentialMaster potentialMaster, Space space,
                                   IRandom random) {
        super(potentialMaster,random,space,1,Double.POSITIVE_INFINITY);//we don't need stepsize-> put 1
        ((MCMoveStepTracker)getTracker()).setTunable(false);
        vCO = space.makeVector();
        vOH = space.makeVector();
    }

    //note that total energy is calculated
    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
        molecule = moleculeSource.getMolecule();

        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();//5

        IAtom c = childList.get(1);
        IAtom sBO = childList.get(3);
        IAtom h = childList.get(4);
        vCO.Ev1Mv2(sBO.getPosition(), c.getPosition());//vector CO
        vOH.Ev1Mv2(h.getPosition(), sBO.getPosition());//vector OH
        double lengthdr13 = vCO.squared();
        
        Vector project = space.makeVector();
        Vector secondaryDirection = space.makeVector();
        project.E(vCO);
        project.TE(vCO.dot(vOH)/lengthdr13);
        secondaryDirection.Ev1Mv2(project,vOH);
        secondaryDirection.TE(2);
        h.getPosition().PE(secondaryDirection);
        for (int k=0; k<numChildren; k++) {
            // shift the whole molecule so that the center of mass (or whatever
            // the position definition uses) doesn't change
        	//newCenter is needed to be changed to oldCenter
            IAtom atomk = childList.get(k);
            atomk.getPosition().PEa1Tv1(-0.2, secondaryDirection);
        }
        uNew = energyMeter.getDataAsScalar();
        return true;
    }

    public void rejectNotify() {
        IAtomList childList = molecule.getChildList();
        int numChildren = childList.size();//5
        IAtom c = childList.get(1);
        IAtom sBO = childList.get(3);
        IAtom h = childList.get(4);
        vCO.Ev1Mv2(sBO.getPosition(), c.getPosition());//vector CO
        vOH.Ev1Mv2(h.getPosition(), sBO.getPosition());//vector OH
        double lengthdr13 = vCO.squared();
        
        Vector project = space.makeVector();
        Vector secondaryDirection = space.makeVector();
        project.E(vCO);
        project.TE(vCO.dot(vOH)/lengthdr13);
        secondaryDirection.Ev1Mv2(project,vOH);
        secondaryDirection.TE(2);
        h.getPosition().PE(secondaryDirection);
        for (int k=0; k<numChildren; k++) {
            IAtom atomk = childList.get(k);
            atomk.getPosition().PEa1Tv1(-0.2, secondaryDirection);
        }
    }

    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld) / temperature);
    }

	
    private static final long serialVersionUID = 1L;
    protected final Vector vCO,vOH;
    
}
