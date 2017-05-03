/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.freeenergy.npath;

import etomica.api.*;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSetSinglet;
import etomica.atom.AtomSource;
import etomica.atom.AtomSourceRandomLeaf;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.atom.iterator.AtomIteratorAtomDependent;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveInsertDeleteLatticeVacancy;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.list.PotentialMasterList;
import etomica.space.ISpace;

/**
 * Created by andrew on 4/11/17.
 */
public class MCMoveAtomSwap extends MCMoveBox {
    protected final AtomIteratorArrayListSimple affectedAtomIterator;
    protected IAtom atom, atom2;
    protected double uOld;
    protected double uNew = Double.NaN;
    protected AtomSource atomSource;
    protected final IRandom random;
    protected ISpace space;
    protected IAtomList atoms;
    protected AtomIteratorAtomDependent atomIterator;
    protected double nbrDistance;
    protected final AtomArrayList nbrList;
    protected final IVectorMutable dr;
    protected final P1ImageHarmonic p1;
    protected final AtomSetSinglet singlet;

    public MCMoveAtomSwap(IRandom random, IPotentialMaster potentialMaster, ISpace _space, P1ImageHarmonic p1) {
        super(potentialMaster);
        this.random = random;
        this.space = _space;
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        perParticleFrequency = true;
        affectedAtomIterator = new AtomIteratorArrayListSimple();
        nbrList = new AtomArrayList();
        dr = space.makeVector();
        this.p1 = p1;
        singlet = new AtomSetSinglet();
    }

    public void setNbrDistance(double newNbrDistance) {
        nbrDistance = newNbrDistance;
    }

    public double getNbrDistance() {
        return nbrDistance;
    }

    /**
     * Method to perform trial move.
     */
    public boolean doTrial() {
        atom = atomSource.getAtom();
        if (atom == null) return false;

        atomIterator.setAtom(atom);
        ((AtomsetIteratorDirectable)atomIterator).setDirection(null);
        atomIterator.reset();
        nbrList.clear();
        IVector pi = atom.getPosition();
        for (IAtom jAtom = ((AtomIterator)atomIterator).nextAtom(); jAtom != null ;jAtom = ((AtomIterator)atomIterator).nextAtom()) {
            dr.Ev1Mv2(pi, jAtom.getPosition());
            box.getBoundary().nearestImage(dr);
            double r2 = dr.squared();
            if (r2 > nbrDistance*nbrDistance) continue;
            nbrList.add(jAtom);
        }
        if (nbrList.getAtomCount() == 0) return false;
        int r = random.nextInt(nbrList.getAtomCount());
        atom2 = nbrList.getAtom(r);

        singlet.atom = atom;
        uOld = 2*p1.energy(singlet);
        singlet.atom = atom2;
        uOld += 2*p1.energy(singlet);

        dr.E(atom.getPosition());
        atom.getPosition().E(atom2.getPosition());
        atom2.getPosition().E(dr);
        return true;
    }//end of doTrial

    /**
     * Returns log of the ratio of the trial probabilities, ln(Tij/Tji) for the
     * states encountered before (i) and after (j) the most recent call to doTrial().
     * Tij is the probability that this move would generate state j from state i, and
     * Tji is the probability that a subsequent call to doTrial would return to state i
     * from state j.
     */
    public double getA() {return 1.0;}

    /**
     * Returns the log of the limiting-distribution probabilities of states, ln(Pj/Pi),
     * for the states encountered before (i) and after (j) the most recent call to
     * doTrial.
     */
    public double getB() {
        singlet.atom = atom;
        uNew = 2*p1.energy(singlet);
        singlet.atom = atom2;
        uNew += 2*p1.energy(singlet);
        return -(uNew - uOld);
    }

    public double energyChange() {return uNew - uOld;}

    /**
     * Method called by IntegratorMC in the event that the most recent trial is accepted.
     */
    public void acceptNotify() {  /* do nothing */
    }

    /**
     * Method called by IntegratorMC in the event that the most recent trial move is
     * rejected.  This method should cause the system to be restored to the condition
     * before the most recent call to doTrial.
     */
    public void rejectNotify() {
        dr.E(atom.getPosition());
        atom.getPosition().E(atom2.getPosition());
        atom2.getPosition().E(dr);
    }

    public AtomIterator affectedAtoms() {
        AtomArrayList list = (AtomArrayList)affectedAtomIterator.getList();
        list.clear();
        list.add(atom);
        list.add(atom2);
        return affectedAtomIterator;
    }

    public void setBox(IBox p) {
        super.setBox(p);
        atomSource.setBox(p);
        atoms = p.getLeafList();

        if (potential instanceof PotentialMasterList) {
            atomIterator = new MCMoveInsertDeleteLatticeVacancy.AtomIteratorNbr(((PotentialMasterList)potential).getNeighborManager(box));
        }
        else if (potential instanceof PotentialMasterCell) {
            atomIterator = new MCMoveInsertDeleteLatticeVacancy.AtomIteratorNbrCell(nbrDistance, ((PotentialMasterCell)potential).getCellAgentManager(), box);
        }
        else {
            // brute force!
            atomIterator = new MCMoveInsertDeleteLatticeVacancy.AtomIteratorBruteForce(box);
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
}
