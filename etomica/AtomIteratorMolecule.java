/*
 * History
 * Created on Dec 30, 2004 by kofke
 */
package etomica;

import etomica.IteratorDirective.Direction;

/**
 * Iterator for the molecules of a single species in a phase.  Can be targeted to
 * give the molecule containing a specific atom (if consistent with the species),
 * or all molecules of the species in a phase.<br>
 * This class is used by PotentialMaster to iterate over molecules for single-body
 * potentials.
 */
public class AtomIteratorMolecule extends AtomIteratorAdapter implements
        AtomsetIteratorMolecule {

    /**
     * @param species species for which molecules are returned as iterates. Only
     * species[0] is relevant, and must not be null.
     */
    public AtomIteratorMolecule(Species[] species) {
        super(new AtomIteratorListSimple());
        listIterator = (AtomIteratorListSimple)iterator;
        this.species = species[0];
        index = this.species.getIndex();
    }

    /**
     * Sets the phase containing the molecules for iteration.
     */
    public void setPhase(Phase phase) {
        speciesAgentNode = (AtomTreeNodeGroup)phase.getAgent(species).node;
        setList();
    }

    /**
     * Specifies molecule returned by this iterator, as the one containing
     * the given target atom.  Only the first element of the array is relevant.
     * If argument is null, of zero length, or if targetAtoms[0] is null, then
     * no target atom is specified and all molecules of the species in the
     * phase will be given on iteration.
     */
    public void setTarget(Atom[] targetAtoms) {
        if(targetAtoms == null || targetAtoms.length == 0) targetAtom = null;
        else targetAtom = targetAtoms[0];
        setList();
    }

    /** 
     * Has no effect, but is included as part of the AtomsetIteratorMolecule interface.
     */
    public void setDirection(Direction direction) {
        //ignore
    }

    /**
     * Configures the list iterator with a list appropriate to the specified
     * phase and target.
     */
    private void setList() {
        if(targetAtom == null) {
            listIterator.setList(speciesAgentNode.childList);
        } else {
            AtomTreeNode moleculeNode = targetAtom.node.childWhereDescendedFrom(speciesAgentNode);
            littleList.clear();
            if(moleculeNode != null) littleList.add(moleculeNode.atom);
            listIterator.setList(littleList);
        }
    }

    private final int index;
    private final AtomIteratorListSimple listIterator;
    private final Species species;
    private final AtomList littleList = new AtomList();
    private AtomTreeNodeGroup speciesAgentNode;
    private Atom targetAtom;
}
