package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.IteratorDirective;
import etomica.Phase;
import etomica.Species;
import etomica.IteratorDirective.Direction;
import etomica.atom.AtomTreeNode;
import etomica.atom.AtomTreeNodeGroup;

/**
 * Gives pairs formed from the molecules of two different species in a phase,
 * taking one molecule of one species with all molecules of the other.
 * Species are specified at construction and cannot be changed afterwards.  The
 * 1-molecule species is identified via the setTarget method, and may be changed
 * from one use of the iterator to the next.
 */

/*
 * History
 * Created on Dec 30, 2004 by kofke
 */

public class ApiInterspecies1A extends AtomsetIteratorAdapter implements
        AtomsetIteratorMolecule {

    /**
     * @param species array of two different, non-null species
     */
    public ApiInterspecies1A(Species[] species) {
        super(new ApiInnerFixed(
                new AtomIteratorSinglet(),
                new AtomIteratorListSimple()));
        if(species[0] == null || species[1] == null) throw new NullPointerException("Constructor of ApiInterspeciesAA requires two non-null species");
        if(species[0] == species[1]) throw new IllegalArgumentException("Constructor of ApiInterspeciesAA requires two different species");
        //arrange so that species0 preceeds species1
        if(species[0].getIndex() < species[1].getIndex()) {
            species0 = species[0];
            species1 = species[1];
        } else {
            species0 = species[1];
            species1 = species[0];
        }
        aiOuter = (AtomIteratorSinglet)((ApiInnerFixed)iterator).getOuterIterator();
        aiInner = (AtomIteratorListSimple)((ApiInnerFixed)iterator).getInnerIterator();
        setPhase(null);
    }

    /** 
     * Configures iterator to return molecules from the set species in the given phase.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        if(phase != null) {
            agentNode0 = (AtomTreeNodeGroup)phase.getAgent(species0).node;
            agentNode1 = (AtomTreeNodeGroup)phase.getAgent(species1).node;
        }
        identifyTargetMolecule();
    }

    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. If the specified direction is consisent with the direction from the
     * target species to the non-target species (as given by their species index --
     * UP is direction from smaller index to larger index) direction, iteration
     * is performed; if specified direction contradicts species direction, no
     * iteration is performed.  Specification of a null direction indicates no
     * limitation, and iteration will be performed if a legitimate target atom 
     * is specified. 
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
        setupIterators();
    }

    /**
     * Sets the target molecule with which all pairs are formed.  Molecule
     * is determined from the first atom of the array, which may be the molecule
     * itself or an atom that is part of it.  If the atom is null or is not 
     * in one of the species given at construction, no iterates will be returned.
     */
    public void setTarget(AtomSet targetAtoms) {
        targetAtom = (Atom)targetAtoms;
        identifyTargetMolecule();
    }

    /**
     * Finds target molecule as indicated by the target atom.  Sets
     * target molecule to null if target atom is null, phase is null, or
     * atom is not part of either species.
     */
    private void identifyTargetMolecule() {
        if(targetAtom == null || phase == null) {
            targetMolecule = null;
        } else {
            AtomTreeNode targetNode = targetAtom.node.childWhereDescendedFrom(agentNode0);

            if(targetNode != null) {    //target is species0
                allowedDirection = IteratorDirective.UP;
                targetMolecule = targetNode.atom();
                aiInner.setList(agentNode1.childList);
                
            } else {                    //target is not species0
                targetNode = targetAtom.node.childWhereDescendedFrom(agentNode1);
                
                if(targetNode != null) {//target is species1
                    allowedDirection = IteratorDirective.DOWN;
                    targetMolecule = targetNode.atom();
                    aiInner.setList(agentNode0.childList);
                    
                } else {                //target not in either species
                    targetMolecule = null;
                }
            }
        }
        setupIterators();
    }
    
    /**
     * Completes setup of iterators, checking that specified direction
     * is consistent with target and species ordering.
     */
    private void setupIterators() {
        if(direction == null || direction == allowedDirection) {
            aiOuter.setAtom(targetMolecule);//targetMolecule may be null here
        } else {
            aiOuter.setAtom(null);
        }
    }
    
    private final AtomIteratorListSimple aiInner;
    private final AtomIteratorSinglet aiOuter;
    private final Species species0, species1;
    
    private AtomTreeNodeGroup agentNode0, agentNode1;
    private Phase phase;
    private Direction direction, allowedDirection;
    private Atom targetAtom, targetMolecule;
}
