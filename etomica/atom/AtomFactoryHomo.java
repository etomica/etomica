package etomica.atom;

import etomica.Atom;
import etomica.Configuration;
import etomica.ConfigurationLinear;
import etomica.Space;

/**
 * Builds an atom group that comprises a set of identically formed atoms or atomgroups.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/23/02 (DAK) changed childFactory from private to protected so that AtomFactoryTree can implement a method to change it
  */
  
public class AtomFactoryHomo extends AtomFactory {
    
    protected AtomFactory childFactory;
    private int atomsPerGroup;
    
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param factory the factory that makes each of the identical children.
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomFactory factory) {
        this(space, sequencerFactory, factory, 1);
    }
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomFactory factory, int atoms) {
        this(space, sequencerFactory, factory, atoms, new ConfigurationLinear(space));
    }
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomFactory factory, 
                            int atoms, Configuration config) {  
        this(space, sequencerFactory, AtomTreeNodeGroup.FACTORY, factory, atoms, config);
    }
 
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param nodeFactory makes nodes for each of the atoms built by this factory
     * @param factory the factory that makes each of the identical children.
     * @param atoms the number of identical children per group (default is 1).
     * @param config the configuration applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomTreeNodeFactory nodeFactory, 
    						AtomFactory factory, int atoms, Configuration config) {
        super(space, sequencerFactory, nodeFactory);
        childFactory = factory;
        atomsPerGroup = atoms;
        configuration = config;
        //set up fields of Group type (can't build sample atoms because factories defined by subclassing this one may not be ready to build at atom at this point)
        
        groupType.childrenAreGroups = factory.isGroupFactory();
    }
    
    public boolean isGroupFactory() {return true;}
    
    /**
     * Constructs a new group.
     */
     public Atom makeAtom() {
        Atom group = newParentAtom();
        for(int i=0; i<atomsPerGroup; i++) {
            Atom childAtom = childFactory.makeAtom();
            childAtom.node.setParent((AtomTreeNodeGroup)group.node);
        }
        configuration.initializeCoordinates(group);
        return group;
     }
     
    /**
     * Returns the subfactory that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory childFactory() {return childFactory;}
        
    /**
     * Specifies the number of child atoms in each atom constructed by this factory.
     * 
     * @param na The new number of atoms per group
     */
    public void setAtomsPerGroup(int na) {
        atomsPerGroup = na;
    }
            
    /**
     * Accessor method for number of child atoms per group constructed.
     */
     public int getAtomsPerGroup() {return atomsPerGroup;}

}//end of AtomFactoryHomo
    
