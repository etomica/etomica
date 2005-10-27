package etomica.atom;

import java.util.Arrays;

import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.data.DataSourceCOM;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactoryNull;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Builds an atom group that comprises a set of differently-formed atoms or
 * atomgroups. Each child atom is constructed by a different atom factory, which
 * are set as an array of atom factories given in the constructor. Position
 * definition is the center-of-mass.
 * 
 * @author David Kofke
 */

public class AtomFactoryHetero extends AtomFactory {

    /**
     * @param factory
     *            array of atom factories, each of which makes a different
     *            child.
     */
    public AtomFactoryHetero(Simulation sim,
            AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType) {
        this(sim.space, sequencerFactory, parentType, new ConformationLinear(sim));
    }

    /**
     * @param factory
     *            the factory that makes each of the identical children.
     * @param atoms
     *            the number of identical children per group (default is 1).
     * @param config
     *            the conformation applied to each group that is built (default
     *            is Linear).
     * @param sequencerFactory
     *            the factory making sequencers used in the groups made by this
     *            factory (default is simple sequencer).
     */
    public AtomFactoryHetero(Space space,
            AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
            Conformation config) {
        this(space, sequencerFactory, parentType, AtomTreeNodeGroup.FACTORY,
                config);
    }

    public AtomFactoryHetero(Space space,
            AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
            AtomTreeNodeFactory nodeFactory, Conformation config) {
        super(new CoordinateFactoryNull(), new AtomTypeGroup(parentType, new DataSourceCOM(space)),
                sequencerFactory, nodeFactory);
        conformation = config;
    }

    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        for (int i = 0; i < childFactory.length; i++) {
            childFactory[i].setSpecies(species);
        }
    }

    /**
     * Constructs a new group.
     */
    public Atom makeAtom() {
        Atom group = newParentAtom();
        AtomTreeNodeGroup node = (AtomTreeNodeGroup) group.node;
        for (int i = 0; i < childFactory.length; i++) {
            for(int j = 0; j < childCount[i]; j++) {
                Atom childAtom = childFactory[i].makeAtom();
                childAtom.node.setParent(node);
            }
        }
        return group;
    }

    /**
     * Sets the factories that make the child atoms of this factory's atom, and
     * specifies the number of child atoms made by each factory.
     * 
     * @param childFactory
     *            array of factories that make the child atoms
     * @param childCount
     *            array giving the number of child atoms each factory makes
     * @throws IllegalArgumentException
     *             if childFactory.length != childCount.length
     */
    public void setChildFactory(AtomFactory[] childFactory, int[] childCount) {
        if (childFactory.length != childCount.length) {
            throw new IllegalArgumentException("childFactory ("
                    + childFactory.length + ") and childCount ("
                    + childCount.length + ") are not the same length");
        }
        this.childFactory = (AtomFactory[]) childFactory.clone();
        this.childCount = (int[]) childCount.clone();
    }

    /**
     * Sets the factories that make the child atoms of this factory's atom,
     * configured so that there is one of each child.
     * 
     * @param childFactory
     */
    public void setChildFactory(AtomFactory[] childFactory) {
        int[] count = new int[childFactory.length];
        Arrays.fill(count, 1);
        setChildFactory(childFactory, count);
    }

    /**
     * 
     * @return
     */

    /**
     * Returns the array of subfactories that produces each of the identical
     * atoms in the group made by this factory.
     */
    public AtomFactory[] getChildFactory() {
        return (AtomFactory[]) childFactory.clone();
    }

    private AtomFactory[] childFactory;
    private int[] childCount;

}//end of AtomFactoryHetero

