package etomica.atom;

import etomica.Atom;
import etomica.Configuration;
import etomica.ConfigurationLinear;
import etomica.Space;
import etomica.atom.AtomTreeNode.Factory;

/**
 * Builds an atom group that comprises a set of differently formed atoms or atomgroups.
 * Each child atom is constructed by a different atom factory, which are set as an
 * array of atomfactories given in the constructor.
 *
 * @author David Kofke
 */
 
 /* History
  * 10/22/02 (DAK) Deleted groupType field which shadowed superclass field, causing errors.
  * 08/26/03 (DAK) Added constructors that permit specification of
  * AtomTreeNode.Factory
  */
  
public class AtomFactoryHetero extends AtomFactory {
    
    private AtomFactory[] childFactory;
    
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
                            Configuration config) {
        this(space, sequencerFactory, AtomTreeNodeGroup.FACTORY, factory, config);
    }
    
	public AtomFactoryHetero(Space space, AtomSequencerFactory sequencerFactory, AtomTreeNodeGroup.Factory nodeFactory,
							AtomFactory[] factory, Configuration config) {
		super(space, sequencerFactory, nodeFactory);
		init(factory, config);
	}
    
    private void init(AtomFactory[] factory, Configuration config) {
        childFactory = factory;
        configuration = config;
        //set up fields of Group type
        for(int i=0; i<factory.length; i++) {
            groupType.childrenAreGroups = factory[i].isGroupFactory();
            if(groupType.childrenAreGroups) break;
        }
    }
    
    public boolean isGroupFactory() {return true;}
    
    /**
     * Constructs a new group.
     */
    public Atom build(Atom group) {
        for(int i=0; i<childFactory.length; i++) {
            childFactory[i].build((AtomTreeNodeGroup)group.node);//builds child atom with group as parent
        }
        configuration.initializeCoordinates(group);
        return group;
    }
    
    /**
     * Returns the array of subfactories that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory[] childFactory() {return childFactory;}
    
    public boolean vetoAddition(Atom a) {return true;} 
        
/*    public void renew(Atom a) {//need an exception in the case a is unrenewable
        if(a.type != groupType) return;  //throw exception
        configuration.initializeCoordinates((AtomGroup)a);
    }       
*/        
}//end of AtomFactoryHomo
    
