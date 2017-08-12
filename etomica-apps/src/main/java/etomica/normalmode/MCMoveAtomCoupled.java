/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move.  Two atoms are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomCoupled extends MCMoveBoxStep {
    
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final MeterPotentialEnergy energyMeter;
    protected final Vector translationVector;
    protected IAtom atom0, atom1;
    protected double uOld;
    protected double uNew;
    protected AtomSource atomSource;
    protected final IRandom random;
    protected IPotentialAtomic[] pairPotential;
    protected final AtomPair pair;
    protected final AtomSetSinglet atomSinglet;
    protected boolean doExcludeNonNeighbors;
    protected int pairPotentialIndex;
    protected IPotentialAtomic constraintPotential;

    public MCMoveAtomCoupled(PotentialMaster potentialMaster, MeterPotentialEnergy energyMeter,
                             IRandom random, Space _space) {
        super(potentialMaster);
        this.random = random;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        this.energyMeter = energyMeter;
        translationVector = _space.makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        energyMeter.setIncludeLrc(false);
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        pair = new AtomPair();
        atomSinglet = new AtomSetSinglet();
        
        setPotential(new IPotentialAtomic[0]);
    }
    
    public void setPotential(IPotentialAtomic[] newPotential) {
    	for(int i=0;i<newPotential.length;i++){
    		if (newPotential[i].nBody() != 2) {
    			throw new RuntimeException("must be a 2-body potential");
    		}
    	}
        pairPotential = newPotential;
    }

    public void setPotential(IPotentialAtomic newPotential) {
    	setPotential(new IPotentialAtomic[] {newPotential});
    }

    public void setConstraint(IPotentialAtomic newConstraintPotential) {
        if (newConstraintPotential.nBody() != 1) {
            throw new RuntimeException("must be a 1-body potential");
        }
        constraintPotential = newConstraintPotential;
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom0 = atomSource.getAtom();
        if (atom0 == null) return false;
        atom1 = atomSource.getAtom();
        if (atom1 == null || atom0 == atom1) return false;
        energyMeter.setTarget(atom0);
        uOld = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom1);
        uOld += energyMeter.getDataAsScalar();
        if(uOld > 1e10) {
            throw new ConfigurationOverlapException(box);
        }
        pair.atom0 = atom0;
        pair.atom1 = atom1;
        pairPotentialIndex = -1;
        if (pairPotential.length > 0 && doExcludeNonNeighbors) {
            IAtomList[] list0 = ((PotentialMasterList)potential).getNeighborManager(box).getDownList(atom0);
            for (int i=0; i<list0.length; i++) {
                if (((AtomArrayList)list0[i]).indexOf(atom1)>-1) {
                	pairPotentialIndex = i;
                    break;
                }
            }
            if (pairPotentialIndex == -1) {
                list0 = ((PotentialMasterList)potential).getNeighborManager(box).getUpList(atom0);
                for (int i=0; i<list0.length; i++) {
                    if (((AtomArrayList)list0[i]).indexOf(atom1)>-1) {
                    	pairPotentialIndex = i;
                        break;
                    }
                }
            }
        }
        if (pairPotentialIndex > -1) {
            uOld -= pairPotential[pairPotentialIndex].energy(pair);
        }

        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom0.getPosition().PE(translationVector);
        atom1.getPosition().ME(translationVector);

        return true;
    }

    public double getChi(double temperature) {
        uNew = energyMeter.getDataAsScalar();
        energyMeter.setTarget(atom0);
        uNew += energyMeter.getDataAsScalar();
        if(!Double.isInfinite(uNew) && pairPotentialIndex > -1){
            uNew -= pairPotential[pairPotentialIndex].energy(pair);
        }
        if (constraintPotential != null) {
            // we could be double-counting the "energy" here (if the atoms are
            // neighbors), but the energy can only be 0 or +inf, so it's OK.
            atomSinglet.atom = atom0;
            constraintPotential.setBox(box);
            uNew += constraintPotential.energy(atomSinglet);
            atomSinglet.atom = atom1;
            uNew += constraintPotential.energy(atomSinglet);
        }

        return Math.exp(-(uNew - uOld) / temperature);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.
     */
    public void rejectNotify() {
        atom0.getPosition().ME(translationVector);
        atom1.getPosition().PE(translationVector);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom0);
        affectedAtomList.add(atom1);
        return affectedAtomIterator;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        atomSource.setBox(p);
        for(int i=0; i< pairPotential.length; i++){
            pairPotential[i].setBox(p);
        }
    }
    
    /**
     * @return Returns the atomSource.
     */
    public AtomSource getAtomSource() {
        return atomSource;
    }
    /**
     * @param source The atomSource to set.
     */
    public void setAtomSource(AtomSource source) {
        atomSource = source;
    }

    /**
     * Configures the move to not explicitly calculate the potential for atom
     * pairs that are not neighbors (as determined by the potentialMaster).
     * 
     * Setting this has an effect only if the potentialMaster is an instance of
     * PotentialMasterList
     */
    public void setDoExcludeNonNeighbors(boolean newDoExcludeNonNeighbors) {
        if (potential == null || !(potential instanceof PotentialMasterList)) {
            throw new RuntimeException("I need a PotentialMasterList to exclude non-neighbors");
        }
        doExcludeNonNeighbors = newDoExcludeNonNeighbors;
    }

    /**
     * Returns true if the move does not explicitly calculate the potential for
     * atom pairs that are not neighbors (as determined by the potentialMaster).
     */
    public boolean getDoExcludeNonNeighbors() {
        return doExcludeNonNeighbors;
    }
}
