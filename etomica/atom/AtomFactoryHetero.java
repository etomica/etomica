package etomica.atom;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.Conformation;
import etomica.ConfigurationLinear;
import etomica.Space;
import etomica.Species;

/**
 * Builds an atom group that comprises a set of differently-formed atoms or atomgroups.
 * Each child atom is constructed by a different atom factory, which are set as an
 * array of atom factories given in the constructor.
 *
 * @author David Kofke
 */
 
public class AtomFactoryHetero extends AtomFactory {
    
    /**
     * @param factory array of atom factories, each of which makes a different child.
     */
	public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomFactory[] factory) {
		this(space, sequencerFactory, factory, new ConfigurationLinear(space));
	}
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     * @param sequencerFactory the factory making sequencers used in the groups made by this factory (default is simple sequencer).
     */
    public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomFactory[] factory, 
                            Conformation config) {
        this(space, sequencerFactory, AtomTreeNodeGroup.FACTORY, factory, config);
    }
    
	public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomTreeNodeFactory nodeFactory,
							AtomFactory[] factory, Conformation config) {
		super(space, new AtomTypeGroup(), sequencerFactory, nodeFactory);
        childFactory = (AtomFactory[])factory.clone();
        configuration = config;
    }
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        for(int i=0; i<childFactory.length; i++) {
            childFactory[i].setSpecies(species);
        }
    }
    /**
     * Constructs a new group.
     */
    public Atom makeAtom() {
        Atom group = newParentAtom();
        AtomTreeNodeGroup node = (AtomTreeNodeGroup)group.node;
        for(int i=0; i<childFactory.length; i++) {
            Atom childAtom = childFactory[i].makeAtom();
            childAtom.node.setParent(node);
        }
        configuration.initializePositions(node.childList);
        return group;
    }
    
    /**
     * Returns the array of subfactories that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory[] childFactory() {return childFactory;}

    private final AtomFactory[] childFactory;

}//end of AtomFactoryHetero
    
