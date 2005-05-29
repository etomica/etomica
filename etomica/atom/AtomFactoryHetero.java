package etomica.atom;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.AtomTreeNodeFactory;
import etomica.AtomTreeNodeGroup;
import etomica.AtomTypeGroup;
import etomica.Conformation;
import etomica.ConformationLinear;
import etomica.Space;
import etomica.Species;
import etomica.data.DataSourceCOM;

/**
 * Builds an atom group that comprises a set of differently-formed atoms or atomgroups.
 * Each child atom is constructed by a different atom factory, which are set as an
 * array of atom factories given in the constructor.  Position definition is the
 * center-of-mass.
 *
 * @author David Kofke
 */
 
public class AtomFactoryHetero extends AtomFactory {
    
    /**
     * @param factory array of atom factories, each of which makes a different child.
     */
	public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType) {
		this(space, sequencerFactory, parentType, new ConformationLinear(space));
	}
    /**
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the conformation applied to each group that is built (default is Linear).
     * @param sequencerFactory the factory making sequencers used in the groups made by this factory (default is simple sequencer).
     */
    public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
                            Conformation config) {
        this(space, sequencerFactory, parentType, AtomTreeNodeGroup.FACTORY, config);
    }
    
	public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
                            AtomTreeNodeFactory nodeFactory, Conformation config) {
		super(space, new AtomTypeGroup(parentType, new DataSourceCOM(space)), sequencerFactory, nodeFactory);
        conformation = config;
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
        return group;
    }
    
    /**
     * @param childFactory The childFactory to set.
     */
    public void setChildFactory(AtomFactory[] childFactory) {
        this.childFactory = (AtomFactory[])childFactory.clone();
    }
    /**
     * Returns the array of subfactories that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory[] getChildFactory() {return (AtomFactory[])childFactory.clone();}

    private AtomFactory[] childFactory;

}//end of AtomFactoryHetero
    
