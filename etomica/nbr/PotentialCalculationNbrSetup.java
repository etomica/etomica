/*
 * History
 * Created on Sep 10, 2004 by kofke
 */
package etomica.nbr;

import etomica.*;

/**
 * PotentialCalculation that performs setup of neighbor lists.  Takes all pair iterates
 * it receives and adds each one of the pair to the other's list of neighbors (which
 * are held in the atoms' sequencers).  This action should be invoked by passing an
 * instance of this class to the PotentialMaster's calculate method (the PotentialMaster
 * should be an instance of PotentialMasterNbr), with an iterator directive that specifies
 * no target atom.
 */
public class PotentialCalculationNbrSetup extends PotentialCalculation {

	public PotentialCalculationNbrSetup() {
		super();
	}

	/**
     * Takes all pair iterates given by the iterator and puts each one of the
     * pair in the other's neighbor list. The second atom in each pair given by
     * the iterator should be uplist of the first atom. Performs no action if
     * the order of the iterator is not exactly 2. The given potential should be
     * a concrete potential, not a potential group (this method is called by the
     * 3-argument PotentialCalculation.doCalculation method after screening out
     * PotentialGroup instances, so this requirement should be met
     * automatically).
     */
	protected void doCalculation(AtomsetIterator iterator, Potential potential) {
		if(iterator.nBody() != 2) return;
		((AtomsetIteratorDirectable)iterator).setDirection(IteratorDirective.UP);
		iterator.reset();
		while(iterator.hasNext()) {
			Atom[] atoms = iterator.next();
			((AtomSequencerNbr)atoms[0].seq).addUpNbr(atoms[1], potential);
			((AtomSequencerNbr)atoms[1].seq).addDownNbr(atoms[0], potential);
		}
	}
}
