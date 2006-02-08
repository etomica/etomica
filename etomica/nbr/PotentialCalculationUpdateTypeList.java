/*
 * History
 * Created on Sep 10, 2004 by kofke
 */
package etomica.nbr;

import java.util.LinkedList;

import etomica.atom.Atom;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.iterator.AtomsetIterator;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.atom.iterator.IteratorDirective;
import etomica.potential.Potential;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialGroup;

/**
 * PotentialCalculation that adds concrete potentials to the NeighborManagerAgents of
 * the AtomTypes to which they apply.  This action should be invoked by passing an
 * instance of this class to the PotentialMaster's calculate method (the PotentialMaster
 * should be an instance of PotentialMasterNbr or PotentialMasterCell).
 */
public class PotentialCalculationUpdateTypeList extends PotentialCalculation {

    public void doCalculation(AtomsetIterator iterator, IteratorDirective id, Potential potential) {    
        if(potential instanceof PotentialGroup) {
            iterator.reset();
            if (iterator.hasNext()) {
                AtomsetIteratorSinglet pretendIterator = new AtomsetIteratorSinglet(iterator.next());
                ((PotentialGroup)potential).calculate(pretendIterator, id, this);
            }
        } else {
            doCalculation(iterator, potential);
        }
    }

    
	/**
     * Takes the first pair given by the iterator and adds the potentials
     * to the atoms' NeighborManagerAgents. Performs no action if
     * the order of the iterator is not 1 or 2. The given potential should be
     * a concrete potential, not a potential group (this method is called by the
     * 3-argument PotentialCalculation.doCalculation method after screening out
     * PotentialGroup instances, so this requirement should be met
     * automatically).
     */
	protected void doCalculation(AtomsetIterator iterator, Potential potential) {
        switch(iterator.nBody()) {
        case 1:
            iterator.reset();
            if(iterator.hasNext()) {
                Atom atom = (Atom)iterator.next();
                potentialAtomTypeList.add(new PotentialAtomTypeWrapper(potential,new AtomType[]{atom.type},iterator));
            }
            break;
        case 2:
            if (iterator instanceof AtomsetIteratorDirectable) {
                ((AtomsetIteratorDirectable)iterator).setDirection(IteratorDirective.Direction.UP);
            }
            iterator.reset();
            if (iterator.hasNext()) {
                AtomPair atoms = (AtomPair)iterator.next();
                potentialAtomTypeList.add(new PotentialAtomTypeWrapper(potential,new AtomType[]{atoms.atom0.type,atoms.atom1.type},iterator));
            }
            break;
        }
	}

    public LinkedList getPotentialsTypeList() {
        return potentialAtomTypeList;
    }

    protected LinkedList potentialAtomTypeList = new LinkedList();

    public static class PotentialAtomTypeWrapper {
        public PotentialAtomTypeWrapper(Potential p, AtomType[] types, AtomsetIterator iterator) {
            potential = p;
            atomTypes = types;
            this.iterator = iterator;
        }
        public final Potential potential;
        public final AtomType[] atomTypes;
        public final AtomsetIterator iterator;
    }
}
