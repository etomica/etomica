/*
 * History
 * Created on Nov 19, 2004 by kofke
 */
package etomica.nbr.cell;

import etomica.Atom;
import etomica.AtomIteratorAtomDependent;
import etomica.AtomIteratorList;
import etomica.AtomIteratorListSimple;
import etomica.AtomIteratorPhaseDependent;
import etomica.AtomLinker;
import etomica.Phase;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.IteratorFactoryCell.CellSequencer;
import etomica.action.AtomsetAction;
import etomica.action.AtomsetCount;
import etomica.action.AtomsetDetect;
import etomica.lattice.AbstractCell;
import etomica.lattice.BravaisLattice;
import etomica.utility.java2.HashMap;

/**
 * @author kofke
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class AtomIteratorNbrCell implements AtomIteratorAtomDependent,
        AtomIteratorPhaseDependent {

    /**
     * 
     */
    public AtomIteratorNbrCell(Species species) {
        this.species = species;
        tabIndex = species.getIndex();
    }

    public void setAtom(Atom atom) {
        this.atom = atom;
        if(atom == null) throw new IllegalArgumentException("Cannot setAtom without referencing an atom");

        AbstractCell referenceCell = ((CellSequencer)atom.seq).cell();
        simpleIterator.setList(referenceCell.neighborManager().neighbors());
    }

    public Atom nextAtom() {
        Atom atom = listIterator.nextAtom();
        if(!listIterator.hasNext()) advanceCell();
        return atom;
    }

    public boolean contains(Atom[] atom) {
        if(atom == null || atom.length == 0) return false;

        AtomsetDetect detector = new AtomsetDetect(atom[0]);
        allAtoms(detector);
        return detector.detectedAtom();
    }

    public boolean hasNext() {
        return listIterator.hasNext();
    }

    public void reset() {
        simpleIterator.reset();//reset cell iterator
        advanceCell();//reset list iterator
    }

    public void unset() {
        simpleIterator.unset();
        listIterator.unset();
    }

    // Moves to next cell that has an iterate
    private void advanceCell() {
        do {
            if(simpleIterator.hasNext()) {//cell iterator
                Atom cell = simpleIterator.nextAtom();
                AtomLinker.Tab[] tabs = (AtomLinker.Tab[])cell.agents[0];
                listIterator.setFirst(tabs[tabIndex]);
            } else {//no more cells at all
                break;
            }
        } while(!listIterator.hasNext());
    }


    public Atom[] next() {
        Atom atom = listIterator.nextAtom();
        if(!listIterator.hasNext()) advanceCell();
        atoms[0] = atom;
        return atoms;
    }

    public Atom[] peek() {
        return listIterator.peek();
    }

    public void allAtoms(AtomsetAction action) {
        if(atom == null) return;
        AbstractCell referenceCell = ((CellSequencer)atom.seq).cell();//cell in which reference atom resides
        HashMap hash = (HashMap)lattice.agents[0];
        int tabIndex = ((Integer)hash.get(speciesAgent)).intValue();//index of tabs in each cell for the basis
        AtomLinker.Tab header = referenceCell.neighborManager().neighbors().header;
        for(AtomLinker e=header.next; e!=header; e=e.next) {//loop over cells
            if(e.atom == null) continue;
            AtomLinker.Tab[] tabs = (AtomLinker.Tab[])e.atom.agents[0];
            AtomLinker next = tabs[tabIndex].next;
            while(next.atom != null) {//loop over atoms in cell
                atoms[0] = next.atom;
                action.actionPerformed(atoms);
                next = next.next;
            }//end while
        }
    }

    public int size() {
        AtomsetCount counter = new AtomsetCount();
        allAtoms(counter);
        return counter.callCount();
    }

    public int nBody() {
        return 1;
    }

    public void setPhase(Phase phase) {
        lattice = getLattice(phase);
        speciesAgent = phase.getAgent(species);
    }

    private int tabIndex;
    private BravaisLattice lattice;
    private Atom atom;
    private Atom[] atoms;
    private Species species;
    private SpeciesAgent speciesAgent;
    private final AtomIteratorList listIterator = new AtomIteratorList();
    private final AtomIteratorListSimple simpleIterator = new AtomIteratorListSimple();

}
