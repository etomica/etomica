/*
 * History
 * Created on Sep 10, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomsetIterator;
import etomica.IteratorDirective;
import etomica.Potential;
import etomica.PotentialGroup;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.atom.iterator.AtomsetIteratorSinglet;
import etomica.potential.PotentialCalculation;

/**
 * PotentialCalculation that adds concrete potentials to the NeighborManagerAgents of
 * the AtomTypes to which they apply.  This action should be invoked by passing an
 * instance of this class to the PotentialMaster's calculate method (the PotentialMaster
 * should be an instance of PotentialMasterNbr or PotentialMasterCell), with an 
 * iterator directive that specifies no target atom.
 */
public class PotentialCalculationAgents extends PotentialCalculation {

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
     * the order of the iterator is not exactly 2. The given potential should be
     * a concrete potential, not a potential group (this method is called by the
     * 3-argument PotentialCalculation.doCalculation method after screening out
     * PotentialGroup instances, so this requirement should be met
     * automatically).
     */
	protected void doCalculation(AtomsetIterator iterator, Potential potential) {
		if(iterator.nBody() != 2) return;
        if (iterator instanceof AtomsetIteratorDirectable) {
            ((AtomsetIteratorDirectable)iterator).setDirection(IteratorDirective.UP);
        }
		iterator.reset();
        if (iterator.hasNext()) {
            AtomPair atoms = (AtomPair)iterator.next();
            Atom atom0 = atoms.atom0;
            Atom atom1 = atoms.atom1;
            atom0.type.getNbrManagerAgent().addPotential(potential);
            if(atom1.type != atom0.type) atom1.type.getNbrManagerAgent().addPotential(potential);
        }
	}
    
    
}
