/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.integrator.mcmove;

import etomica.api.*;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.util.Debug;

/**
 * An MC Move for cluster simulations that "wiggles" for acetic acid.
 * 
 * @author Hye Min Kim
 */
public class MCMoveWiggleAceticAcid extends MCMoveMolecule {

    public MCMoveWiggleAceticAcid(Simulation sim, PotentialMaster potentialMaster, Space _space) {
    	this(potentialMaster,sim.getRandom(), 0.1, _space);
    }
    
    public MCMoveWiggleAceticAcid(PotentialMaster potentialMaster,
            IRandom random, double stepSize, Space _space) {
        super(potentialMaster,random,_space, stepSize,Double.POSITIVE_INFINITY);
        this.space = _space;
        bondedAtoms = new int[]{1,-1,1,-1,3};//0(Ch3) and 2(dBO) are bonded to 1 (C), 4(H) is bonded to 3 (O)
        setStepSizeMax(Math.PI);
        work1 = _space.makeVector();
        work2 = _space.makeVector();
        work3 = _space.makeVector();
        translationVectors = space.makeVector();
    }

    //note that total energy is calculated
    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
        molecule = moleculeSource.getMolecule();

        IAtomList childList = molecule.getChildList();
        int numChildren = childList.getAtomCount();//5

        int j = random.nextInt(3)*2;
        selectedAtoms = j;//0,2,4
        IAtom selectedAtom = childList.getAtom(j);
        Vector position = selectedAtom.getPosition();
        translationVectors.Ea1Tv1(-1,position);
        double oldBondLength1 = 0, oldBondLength2 = 0;

        // this puts atom j in a random orientation without changing the bond length

        //work1 is the current vector from the bonded atom to atom j
        work1.E(position);
        work1.ME(childList.getAtom(bondedAtoms[j]).getPosition());
        position.E(childList.getAtom(bondedAtoms[j]).getPosition());
        double bondLength = Math.sqrt(work1.squared());
        if (Debug.ON && Debug.DEBUG_NOW) {
            oldBondLength1 = bondLength;
        }
        //work2 is a vector perpendicular to work1.  it can be any 
        //perpendicular vector, but that just makes it harder!
        if (work1.getX(0)*work1.getX(0) < 0.5*bondLength*bondLength) {
            // if work1 doesn't point in the X direction (mostly) then
            // find a vector in the plane containing the X axis and work1
            double a = -work1.getX(0)/bondLength;
            work2.Ea1Tv1(a,work1);
            work2.setX(0,work2.getX(0)+bondLength);
        }
        else {
            // work1 does point in the X direction (mostly) so
            // find a vector in the plane containing the Y axis and work1
            double a = -work1.getX(1)/bondLength;
            work2.Ea1Tv1(a,work1);
            work2.setX(1,work2.getX(1)+bondLength);
        }
        //normalize
        work2.TE(bondLength/Math.sqrt(work2.squared()));
        //work3 is a vector normal to both work1 and work2
        work3.E(work1);
        work3.XE(work2);
        work3.TE(bondLength/Math.sqrt(work3.squared()));
        
        double phi = (random.nextDouble()-0.5)*Math.PI;
        work2.TE(Math.cos(phi));
        work2.PEa1Tv1(Math.sin(phi),work3);

        double theta = (random.nextDouble()-0.5)*stepSize;
        position.PEa1Tv1(Math.cos(theta),work1);
        position.PEa1Tv1(Math.sin(theta),work2);

        translationVectors.PE(position);
        work1.E(translationVectors);
        work1.TE(1.0/childList.getAtomCount());
        for (int k=0; k<childList.getAtomCount(); k++) {//handling COM
            childList.getAtom(k).getPosition().ME(work1);
        }
        if (Debug.ON && Debug.DEBUG_NOW) {
            if (j > 0) {
                work1.Ev1Mv2(position, childList.getAtom(j-1).getPosition());
                bondLength = Math.sqrt(work1.squared());
                if (Math.abs(bondLength - oldBondLength1)/oldBondLength1 > 0.000001) {
                    throw new IllegalStateException("wiggle "+" "+j+" bond length should be close to "+oldBondLength1+" ("+bondLength+")");
                }
            }
            if (j < numChildren-1) {
                work1.Ev1Mv2(position, childList.getAtom(j+1).getPosition());
                bondLength = Math.sqrt(work1.squared());
                double oldBondLength = oldBondLength2 == 0 ? oldBondLength1 : oldBondLength2;
                if (Math.abs(bondLength - oldBondLength)/oldBondLength > 0.000001) {
                    throw new IllegalStateException("wiggle "+" "+j+" bond length should be close to "+oldBondLength+" ("+bondLength+")");
                }
            }
        }

        uNew = energyMeter.getDataAsScalar();
        return true;
    }

    public void rejectNotify() {
        IAtomList childList = molecule.getChildList();
        work1.E(translationVectors);
        work1.TE(1.0/childList.getAtomCount());
        for (int k=0; k<childList.getAtomCount(); k++) {
            childList.getAtom(k).getPosition().PE(work1);//undo COM
        }
        childList.getAtom(selectedAtoms).getPosition().ME(translationVectors);
    }

    public double getB() {
        return -(uNew - uOld);
    }
	
    private static final long serialVersionUID = 1L;
    protected int selectedAtoms;
    protected int[] bondedAtoms;
    protected final Vector work1, work2, work3;
    protected final Vector translationVectors;
    protected final Space space;
}
