package etomica.atom.iterator;

import java.io.Serializable;

import etomica.action.AtomsetAction;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.SpeciesAgent;
import etomica.atom.iterator.IteratorDirective.Direction;
import etomica.phase.Phase;
import etomica.species.Species;

/**
 * Gives pairs formed from the molecules of two different species in a phase,
 * taking one molecule of one species with all molecules of the other. Species
 * are specified at construction and cannot be changed afterwards. The
 * 1-molecule species is identified via the setTarget method, and may be changed
 * from one use of the iterator to the next.
 * <p>
 * Atoms in iterate pairs are ordered such that the first atom is of the first
 * species given in the constructor, and the second atom is of the second
 * species given in the constructor, regardless of the specification of the
 * target.
 * <p>
 * Direction may be specified, and is interpreted as follows:
 * <ul>
 * <li>For null direction, no restriction on iterates.
 * <li>For UP direction, iterates are given if index of target-atom species is
 * less than index of other species.
 * <li>For DOWN direction, iterates are given if index of target-atom species
 * is greater than index of other species.
 * </ul>
 */

public class ApiInterspecies1A implements AtomPairIterator, AtomsetIteratorPDT,
        Serializable {

    /**
     * Sorts given array of species according to species index, then constructs iterator 
     * such that atom0 of the pair iterates is in (sorted) species[0], and atom1 is in 
     * (sorted) species[1], regardless of which is specified via setTarget.  Thus the
     * species index of atom0 is less than that of atom1, for all iterates.
     * 
     * @param species
     *            array of two different, non-null species
     * 
     * @throws IllegalArgumentException
     *             is species array is not of length 2 or if species in array
     *             refer to the same instance
     * @throws NullPointerException
     *             if species array is null or if either species in array is
     *             null
     */
    public ApiInterspecies1A(Species[] species) {
        super();
        if (species.length != 2) {
            throw new IllegalArgumentException(
                    "Constructor of ApiInterspecies1A requires an array of two species");
        }
        if (species[0] == null || species[1] == null) {
            throw new NullPointerException(
                    "Constructor of ApiInterspeciesAA requires two non-null species");
        }
        if (species[0] == species[1]) {
            throw new IllegalArgumentException(
                    "Constructor of ApiInterspeciesAA requires two different species");
        }
        aiOuter = new AtomIteratorSinglet();
        aiInner = new AtomIteratorArrayListSimple();
        apiUp = new ApiInnerFixed(aiOuter, aiInner, false);
        apiDown = new ApiInnerFixed(aiOuter, aiInner, true);
        iterator = apiUp;

        // we need to sort these.  we'll do that once we have the phase
        species0 = species[0];
        species1 = species[1];
        setPhase(null);
    }

    /**
     * Configures iterator to return molecules from the set species in the given
     * phase.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        agent0 = phase.getAgent(species0);
        agent1 = phase.getAgent(species1);
        if (agent0.getIndex() > agent1.getIndex()) {
            // species were out of order.  swap them
            Species tempSpecies = species0;
            species0 = species1;
            species1 = tempSpecies;
            SpeciesAgent tempAgent = agent0;
            agent0 = agent1;
            agent1 = tempAgent;
        }
        identifyTargetMolecule();
    }

    /**
     * Indicates allowed direction for iteration, relative to specified target
     * atom. If the specified direction is consisent with the direction from the
     * target species to the non-target species (as given by their species index --
     * UP is direction from smaller index to larger index) direction, iteration
     * is performed; if specified direction contradicts species direction, no
     * iteration is performed. Specification of a null direction indicates no
     * limitation, and iteration will be performed if a legitimate target atom
     * is specified.
     */
    public void setDirection(Direction direction) {
        this.direction = direction;
        setupIterators();
    }

    /**
     * Sets the target molecule with which all pairs are formed. Molecule is
     * determined from the atom specified by the newTargetAtom (which must not
     * be null), which may be the molecule itself or an atom that is part
     * of it. If the atom is not in one of the species given at
     * construction, no iterates will be returned.
     * 
     * @throws NullPointerException
     *             if targetAtom is null
     */
    public void setTarget(IAtom newTargetAtom) {
        if (newTargetAtom == null) {
            throw new NullPointerException("target atom must not be null");
        }
        targetAtom = newTargetAtom;
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
        }
        else {
            targetMolecule = targetAtom.getChildWhereDescendedFrom(agent0);
            if (targetMolecule != null) {
                //target is species0
                allowedDirection = IteratorDirective.Direction.UP;
                iterator = apiUp;
                aiInner.setList(agent1.getChildList());
            }
            else {
                //target is not species0, try species1
                targetMolecule = targetAtom.getChildWhereDescendedFrom(agent1);

                if (targetMolecule != null) {//target is species1
                    allowedDirection = IteratorDirective.Direction.DOWN;
                    iterator = apiDown;
                    aiInner.setList(agent0.getChildList());
                }
            }
        }
        setupIterators();
    }

    /**
     * Completes setup of iterators, checking that specified direction is
     * consistent with target and species ordering.
     */
    private void setupIterators() {
        if (direction == null || direction == allowedDirection) {
            aiOuter.setAtom(targetMolecule);//targetMolecule may be null here
        } else {
            aiOuter.setAtom(null);
        }
    }

    public AtomPair nextPair() {
        return iterator.nextPair();
    }

    public void allAtoms(AtomsetAction action) {
        iterator.allAtoms(action);
    }
    
    public boolean hasNext() {
        return iterator.hasNext();
    }
    
    public int nBody() {
        return 2;
    }
    
    public AtomSet next() {
        return iterator.nextPair();
    }
    
    public void reset() {
        iterator.reset();
    }

    public int size() {
        return iterator.size();
    }
    
    public void unset() {
        iterator.unset();
    }
    
    private static final long serialVersionUID = 1L;
    private final AtomIteratorArrayListSimple aiInner;
    private final AtomIteratorSinglet aiOuter;
    private Species species0, species1;
    private final ApiInnerFixed apiUp, apiDown;
    private ApiInnerFixed iterator;
    private IteratorDirective.Direction direction, allowedDirection;
    private SpeciesAgent agent0, agent1;
    private Phase phase;
    private IAtom targetAtom, targetMolecule;
}
