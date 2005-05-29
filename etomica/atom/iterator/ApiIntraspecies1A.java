package etomica.atom.iterator;

import etomica.Atom;
import etomica.AtomSet;
import etomica.AtomTreeNode;
import etomica.AtomTreeNodeGroup;
import etomica.Phase;
import etomica.Species;
import etomica.IteratorDirective.Direction;

/**
 * Gives pairs formed from the molecules of a species in a phase,
 * taking one molecule the species with all of its other molecules.
 * Species is specified at construction and cannot be changed afterwards.
 */

/*
 * History
 * Created on Dec 30, 2004 by kofke
 */

public class ApiIntraspecies1A extends AtomPairIteratorAdapter implements
        AtomsetIteratorMolecule {

    /**
     * @param species species whose molecules will form the pair iterates
     */
    public ApiIntraspecies1A(Species species) {
        this(new Species[] {species, species});
    }
    
    /**
     * @param species array of two non-null elements referencing the same species instance
     */
    public ApiIntraspecies1A(Species[] species) {
        super(new ApiInnerVariable(
                new AtomIteratorSinglet(),
                new AtomIteratorSequencerList()));
        if(species == null || species.length < 1 || species[0] == null) throw new NullPointerException("Constructor of ApiIntraspecies1A requires two non-null species references to the same instance");
        if(species[0] != species[1]) throw new IllegalArgumentException("Constructor of ApiIntraspecies1A requires references to the same species instance");
        this.species = species[0];

        aiOuter = (AtomIteratorSinglet)((ApiInnerVariable)iterator).getOuterIterator();
        aiInner = (AtomIteratorSequencerList)((ApiInnerVariable)iterator).getInnerIterator();
        aiInner.setNumToSkip(1);
        setPhase(null);
    }
    /** 
     * Configures iterator to return molecules from the set species in the given phase.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        if(phase != null) {
            agentNode = (AtomTreeNodeGroup)phase.getAgent(species).node;
            identifyTargetMolecule();
        } else {
            targetMolecule = null;
            aiOuter.setAtom(null);
            aiInner.setAtom(null);
        }
    }

    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. Specification of a null direction indicates iteration in both directions
     * relative to the target. 
     */
    public void setDirection(Direction direction) {
        aiInner.setDirection(direction);
    }

    /**
     * Sets the target molecule with which all pairs are formed.  Molecule
     * is determined from the first atom of the array, which may be the molecule
     * itself or an atom that is part of it.  If the atom is null or is not 
     * in one of the species given at construction, no iterates will be returned.
     */
    public void setTarget(AtomSet targetAtoms) {
        if(targetAtoms.count() != 1) throw new IllegalArgumentException("1A iterator must have exactly one target atom");
        targetAtom = targetAtoms.getAtom(0);
        identifyTargetMolecule();
    }

    /**
     * Finds target molecule as indicated by the target atom.  Sets
     * target molecule to null if target atom is null, phase is null, or
     * atom is not part of either species.
     */
    private void identifyTargetMolecule() {
        if(phase == null) {
            targetMolecule = null;
        } else {
            AtomTreeNode targetNode = targetAtom.node.childWhereDescendedFrom(agentNode);
            targetMolecule = (targetNode != null) ? targetNode.atom() : null;
        }
        //targetMolecule may be null here
        aiOuter.setAtom(targetMolecule);
        aiInner.setAtom(targetMolecule);
    }
    
    private final AtomIteratorSequencerList aiInner;
    private final AtomIteratorSinglet aiOuter;
    private final Species species;
    
    private AtomTreeNodeGroup agentNode;
    private Phase phase;
    private Atom targetAtom, targetMolecule;
}
