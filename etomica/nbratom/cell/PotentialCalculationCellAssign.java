/*
 * History
 * Created on Sep 10, 2004 by kofke
 */
package etomica.nbratom.cell;

import etomica.Atom;
import etomica.AtomPair;
import etomica.AtomPairIterator;
import etomica.AtomsetIterator;
import etomica.IteratorDirective;
import etomica.Potential;
import etomica.atom.iterator.AtomsetIteratorDirectable;
import etomica.potential.PotentialCalculation;

/**
 * PotentialCalculation that performs setup of neighbor lists.  Takes all pair iterates
 * it receives and adds each one of the pair to the other's list of neighbors (which
 * are held in the atoms' sequencers).  This action should be invoked by passing an
 * instance of this class to the PotentialMaster's calculate method (the PotentialMaster
 * should be an instance of PotentialMasterNbr), with an iterator directive that specifies
 * no target atom.
 */
public class PotentialCalculationCellAssign extends PotentialCalculation {

	public PotentialCalculationCellAssign(NeighborCellManager cellManager) {
		super();
        this.cellManager = cellManager;
	}

	/**
     * Takes all pair iterates given by the iterator and puts each one of the
     * pair in the other's neighbor list. Performs no action if
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
            AtomPair atoms = (AtomPair)iterator.peek();
            Atom atom0 = atoms.atom0;
            Atom atom1 = atoms.atom1;
            atom0.type.getNbrManagerAgent().addPotential(potential);
            if(atom1.type != atom0.type) atom1.type.getNbrManagerAgent().addPotential(potential);
        }
		while(iterator.hasNext()) {
			AtomPair atoms = ((AtomPairIterator)iterator).nextPair();
            Atom atom0 = atoms.atom0;
            Atom atom1 = atoms.atom1;
            if(((AtomSequencerCell)atom0.seq).cell == null) {
                cellManager.assignCell(atom0);
            }
            if(((AtomSequencerCell)atom1.seq).cell == null) {
                cellManager.assignCell(atom1);
            }
		}
	}
    
    private final NeighborCellManager cellManager;
    
}
