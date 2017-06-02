/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.*;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.Vector3D;
import etomica.util.Debug;

/**
 * An MC Move for cluster simulations that "wiggles" a chain molecule.  If the 
 * first or last atom in the chain is chosen, it is moved to a new position 
 * with the same bond length as before, but perturbed by some angle from its
 * original position.  If an Atom in the middle of the chain is chosen, a 
 * crankshaft move is performed that maintains its distances with its 
 * neighbors.  If a middle Atom has a bond angle too close to 180 degrees
 * (such that rotation does nothing) the Atom is not moved at all.
 * In each doTrial, wiggle moves are attempted on all molecules in the box. 
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterWiggleAceticAcid extends MCMoveMolecule {

    public MCMoveClusterWiggleAceticAcid(Simulation sim, PotentialMaster potentialMaster, Space _space) {
    	this(potentialMaster,sim.getRandom(), 0.1, _space);//0.1 rad wiggle move
    }
    
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in
     * box should be at least one greater than this value (greater
     * because first atom is never moved)
     */
    public MCMoveClusterWiggleAceticAcid(PotentialMaster potentialMaster,
            IRandom random, double stepSize, Space _space) {
        super(potentialMaster,random,_space, stepSize,Double.POSITIVE_INFINITY);
        this.space = _space;
        bondedAtoms = new int[]{1,-1,1,-1,3};//0(Ch3) and 2(dBO) are bonded to 1 (C), 4(H) is bonded to 3 (O)
        setStepSizeMax(Math.PI);
        energyMeter = new MeterPotentialEnergy(potential);
        work1 = _space.makeVector();
        work2 = _space.makeVector();
        work3 = _space.makeVector();
    }

    public void setBox(Box p) {
        super.setBox(p);
        selectedAtoms = new int[box.getMoleculeList().getMoleculeCount()];
        translationVectors = new Vector3D[box.getMoleculeList().getMoleculeCount()];
        for (int i=0; i<translationVectors.length; i++) {
            translationVectors[i] = space.makeVector();
        }
        energyMeter.setBox(p);
    }
    
    public void setSpecies(ISpecies newSpecies) {
        species = newSpecies;
    }

    //note that total energy is calculated
    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<moleculeList.getMoleculeCount(); i++) {
            if (species != null && moleculeList.getMolecule(i).getType() != species) {
                continue;
            }
            IAtomList childList = moleculeList.getMolecule(i).getChildList();
            int numChildren = childList.getAtomCount();

            int j = random.nextInt(3)*2;
            selectedAtoms[i] = j;//0,2,4
            IAtom selectedAtom = childList.getAtom(j);
            Vector position = selectedAtom.getPosition();
            translationVectors[i].Ea1Tv1(-1,position);
            double oldBondLength1 = 0, oldBondLength2 = 0;

            // this puts atom j in a random orientation without changing
            // the bond length

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

            translationVectors[i].PE(position);
            work1.E(translationVectors[i]);
            work1.TE(1.0/childList.getAtomCount());
            for (int k=0; k<childList.getAtomCount(); k++) {//handling COM
                childList.getAtom(k).getPosition().ME(work1);
            }
            if (Debug.ON && Debug.DEBUG_NOW) {
                if (j > 0) {
                    work1.Ev1Mv2(position, childList.getAtom(j-1).getPosition());
                    bondLength = Math.sqrt(work1.squared());
                    if (Math.abs(bondLength - oldBondLength1)/oldBondLength1 > 0.000001) {
                        throw new IllegalStateException("wiggle "+i+" "+j+" bond length should be close to "+oldBondLength1+" ("+bondLength+")");
                    }
                }
                if (j < numChildren-1) {
                    work1.Ev1Mv2(position, childList.getAtom(j+1).getPosition());
                    bondLength = Math.sqrt(work1.squared());
                    double oldBondLength = oldBondLength2 == 0 ? oldBondLength1 : oldBondLength2;
                    if (Math.abs(bondLength - oldBondLength)/oldBondLength > 0.000001) {
                        throw new IllegalStateException("wiggle "+i+" "+j+" bond length should be close to "+oldBondLength+" ("+bondLength+")");
                    }
                }
            }
        }
        ((BoxCluster)box).trialNotify();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i=0; i<selectedAtoms.length; i++) {
            if (species != null && moleculeList.getMolecule(i).getType() != species) continue;
            IAtomList childList = moleculeList.getMolecule(i).getChildList();
            work1.E(translationVectors[i]);
            work1.TE(1.0/childList.getAtomCount());
            for (int k=0; k<childList.getAtomCount(); k++) {
                childList.getAtom(k).getPosition().PE(work1);//undo COM
            }
            childList.getAtom(selectedAtoms[i]).getPosition().ME(translationVectors[i]);
        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getB() {
        return -(uNew - uOld);
    }
    
    public double getA() {
        return wNew/wOld;
    }
	
    private static final long serialVersionUID = 1L;
    protected final MeterPotentialEnergy energyMeter;
    protected int[] selectedAtoms;
    protected int[] bondedAtoms;
    protected final Vector work1, work2, work3;
    protected Vector[] translationVectors;
    protected double wOld, wNew;
    protected final Space space;
    protected ISpecies species;
}
