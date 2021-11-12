/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.IPotentialField;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Standard Monte Carlo atom-displacement trial move.  Two atoms are moved at a
 * time in such a way that the geometric center of the system is not changed.
 *
 * @author Andrew Schultz
 */
public class MCMoveAtomCoupledFasterer extends MCMoveBoxStep {

    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected final AtomArrayList affectedAtomList;
    protected final PotentialCompute potentialCompute;
    protected final Vector translationVector;
    protected IAtom atom1, atom2;
    protected double uOld;
    protected double uNew;
    protected AtomSource atomSource;
    protected final IRandom random;
    protected Vector oldPosition1, oldPosition2;
    protected IPotentialField constraintPotential;
    protected boolean callComputeManyAtoms;

    public MCMoveAtomCoupledFasterer(PotentialCompute potentialCompute, IRandom random, Space _space) {
        super();
        this.potentialCompute = potentialCompute;
        this.random = random;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        translationVector = _space.makeVector();
        setStepSizeMax(0.5);
        setStepSizeMin(0.0);
        setStepSize(0.1);
        perParticleFrequency = true;
        affectedAtomList = new AtomArrayList(2);
        affectedAtomIterator = new AtomIteratorArrayListSimple(affectedAtomList);
        oldPosition1 = _space.makeVector();
        oldPosition2 = _space.makeVector();
    }

    /**
     * Informs the move that it should call computeManyAtoms instead of
     * computeManyAtomsOld.  This will allow systems with a neighbor range -
     * potential range mismatch to work and will also work with a lattice sum.
     * Calling computeManyAtoms may be slower and may break other MC moves that
     * rely on the PotentialCompute keeping track of the energy of every atom.
     */
    public void setCallComputeManyAtoms(boolean doCallComputeManyAtoms) {
        callComputeManyAtoms = doCallComputeManyAtoms;
    }

    public void setConstraint(IPotentialField newConstraintPotential) {
        if (newConstraintPotential.nBody() != 1) {
            throw new RuntimeException("must be a 1-body potential");
        }
        constraintPotential = newConstraintPotential;
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom1 = atomSource.getAtom();
        if (atom1 == null) return false;
        atom2 = atomSource.getAtom();
        if (atom2 == null || atom1 == atom2) return false;
        if (callComputeManyAtoms) {
            uOld = potentialCompute.computeManyAtoms(atom1, atom2);
        }
        else {
            uOld = potentialCompute.computeManyAtomsOld(atom1, atom2);
        }
        if(uOld > 1e10) {
            throw new ConfigurationOverlapException(box);
        }
        oldPosition1.E(atom1.getPosition());
        oldPosition2.E(atom2.getPosition());

        translationVector.setRandomCube(random);
        translationVector.TE(stepSize);
        atom1.getPosition().PE(translationVector);
        atom2.getPosition().ME(translationVector);

        return true;
    }

    public double getChi(double temperature) {
        uNew = potentialCompute.computeManyAtoms(atom1, atom2);

        if (constraintPotential != null) {
            // we could be double-counting the "energy" here (if the atoms are
            // neighbors), but the energy can only be 0 or +inf, so it's OK.
            uNew += constraintPotential.u(atom1);
            uNew += constraintPotential.u(atom2);
        }

        return Math.exp(-(uNew - uOld) / temperature);
    }
    
    public double energyChange() {return uNew - uOld;}
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {
        if (!callComputeManyAtoms) {
            potentialCompute.processAtomU(1);
            // put it back, then compute old contributions to energy
            atom1.getPosition().E(oldPosition1);
            atom2.getPosition().E(oldPosition2);

            potentialCompute.computeManyAtoms(atom1, atom2);

            atom1.getPosition().PE(translationVector);
            atom2.getPosition().ME(translationVector);

            potentialCompute.processAtomU(-1);
        }
    }
    
    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.
     */
    public void rejectNotify() {
        atom1.getPosition().E(oldPosition1);
        atom2.getPosition().E(oldPosition2);
    }
        
    
    public AtomIterator affectedAtoms() {
        affectedAtomList.clear();
        affectedAtomList.add(atom1);
        affectedAtomList.add(atom2);
        return affectedAtomIterator;
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        atomSource.setBox(p);
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
}
