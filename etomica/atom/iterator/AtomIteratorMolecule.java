/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica.atom.iterator;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNode;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Iterator for the molecules of a single species in a phase.  Can be targeted to
 * give the molecule containing a specific atom (if consistent with the species),
 * or all molecules of the species in a phase.<br>
 * This class is used by PotentialMaster to iterate over molecules for single-body
 * potentials.
 */
public class AtomIteratorMolecule extends AtomIteratorAdapter implements
        AtomsetIteratorPDT, AtomIteratorPhaseDependent {

    /**
     * @param species species for which molecules are returned as iterates. Only
     * species[0] is relevant, and must not be null.
     */
    public AtomIteratorMolecule(Species[] species) {
        super(new AtomIteratorArrayListSimple());
        listIterator = (AtomIteratorArrayListSimple)iterator;
        this.species = species[0];
        setList();
    }

    /**
     * Sets the phase containing the molecules for iteration. A null
     * phase conditions iterator to give no iterates.
     */
    public void setPhase(Phase phase) {
        if(phase == null) speciesAgentNode = null;
        else speciesAgentNode = (AtomTreeNodeGroup)phase.getAgent(species).node;
        setList();
    }

    /**
     * Specifies molecule returned by this iterator, as the one containing
     * the given target atom.  Only the first element of the array is relevant.
     * If argument is null, of zero length, or if targetAtoms[0] is null, then
     * no target atom is specified and all molecules of the species in the
     * phase will be given on iteration.
     * 
     * @throws NullPointerException
     *          if targetAtoms is null
     * @throws IllegalArgumentException
     *          if targetAtoms.count() is not 0 or 1
     */
    public void setTarget(AtomSet targetAtoms) {
        switch(targetAtoms.count()) {
        case 0: 
            targetAtom = null;
            break;
        case 1:
            targetAtom = targetAtoms.getAtom(0);
            break;
        default:
            throw new IllegalArgumentException("Can specify at most one target atom to iterator");
        }
        setList();
    }

    /** 
     * Has no effect, but is included as part of the AtomsetIteratorPDT interface.
     */
    public void setDirection(Direction direction) {
        //ignore
    }

    /**
     * Configures the list iterator with a list appropriate to the specified
     * phase and target.
     */
    private void setList() {
        //no phase is specified
        if(speciesAgentNode == null) {
            listIterator.setList(null);
            
        //no target -- iterate all molecules of species
        } else if(targetAtom == null) {
            listIterator.setList(speciesAgentNode.childList);
        
        //target specified -- give it as only iterate if descended from species
        } else {
            AtomTreeNode moleculeNode = targetAtom.node.childWhereDescendedFrom(speciesAgentNode);
            littleList.clear();
            if(moleculeNode != null) littleList.add(moleculeNode.atom());
            listIterator.setList(littleList);
        }
    }

    private final AtomIteratorArrayListSimple listIterator;
    private final Species species;
    private final AtomArrayList littleList = new AtomArrayList();
    private AtomTreeNodeGroup speciesAgentNode;
    private Atom targetAtom;
}
