package etomica.atom.iterator;

import etomica.atom.Atom;
import etomica.atom.AtomSet;
import etomica.atom.AtomTreeNode;
import etomica.atom.AtomTreeNodeGroup;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Gives pairs formed from the molecules of a species in a phase, taking one
 * molecule the species with all of its other molecules. Species is specified at
 * construction and cannot be changed afterwards.
 */

/*
 * History Created on Dec 30, 2004 by kofke
 */

public class ApiIntraspecies1A extends ApiSequence1A implements
        AtomsetIteratorMolecule {

    /**
     * @param species
     *            species whose molecules will form the pair iterates
     */
    public ApiIntraspecies1A(Species species) {
        this(new Species[] { species, species });
    }

    /**
     * @param species
     *            array of two non-null elements referencing the same species
     *            instance
     * @throws NullPointerException
     *             if species or one of its elements is null
     * @throws IllegalArgumentException
     *             if species is not a length-2 array or if its elements refer
     *             to different species instances
     */
    public ApiIntraspecies1A(Species[] species) {
        super();
        if (species.length != 2)
            throw new IllegalArgumentException(
                    "Constructor of ApiIntraspecies1A requires two references to the same species instance");
        if (species == null || species[0] == null || species[1] == null)
            throw new NullPointerException(
                    "Constructor of ApiIntraspecies1A requires two non-null species references to the same instance");
        if (species[0] != species[1])
            throw new IllegalArgumentException(
                    "Constructor of ApiIntraspecies1A requires references to the same species instance");
        this.species = species[0];

        setPhase(null);
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * phase. No iterates are given if phase is null.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        if (phase != null) {
            agentNode = (AtomTreeNodeGroup) phase.getAgent(species).node;
            identifyTargetMolecule();
        } else {
            targetMolecule = null;
            setAtom(null);
        }
    }

    /**
     * Sets the target molecule with which all pairs are formed. Molecule is
     * determined from the atom specified by the atomSet (which must have
     * count() == 1), which may be the molecule itself or an atom that is part
     * of it. If the atom is null or is not in the species given at
     * construction, no iterates will be returned.
     * 
     * @throws NullPointerException
     *             if targetAtoms is null
     * @throws IllegalArgumentException
     *             if targetAtoms.count() is not 1
     */
    public void setTarget(AtomSet targetAtoms) {
        if (targetAtoms.count() != 1)
            throw new IllegalArgumentException(
                    "1A iterator must have exactly one target atom (which may be null); count of given target is "
                            + targetAtoms.count());
        targetAtom = targetAtoms.getAtom(0);
        identifyTargetMolecule();
    }

    /**
     * Finds target molecule as indicated by the target atom. Sets target
     * molecule to null if target atom is null, phase is null, or atom is not
     * part of either species.
     */
    private void identifyTargetMolecule() {
        if (phase == null || targetAtom == null) {
            targetMolecule = null;
        } else {
            AtomTreeNode targetNode = targetAtom.node
                    .childWhereDescendedFrom(agentNode);
            targetMolecule = (targetNode != null) ? targetNode.atom() : null;
        }
        //targetMolecule may be null here
        setAtom(targetMolecule);
    }

    private final Species species;

    private AtomTreeNodeGroup agentNode;
    private Phase phase;
    private Atom targetAtom, targetMolecule;

}
