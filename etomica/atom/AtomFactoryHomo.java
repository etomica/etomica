package etomica.atom;

import etomica.Atom;
import etomica.AtomFactory;
import etomica.Conformation;
import etomica.ConfigurationLinear;
import etomica.Space;
import etomica.Species;

/**
 * Builds an atom group that comprises a set of identically formed atoms or atom groups.
 *
 * @author David Kofke
 */
 
 public class AtomFactoryHomo extends AtomFactory {
    
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
                            int atoms, Conformation config) {  
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
    						AtomFactory factory, int atoms, Conformation config) {
        super(space, new AtomTypeGroup(), sequencerFactory, nodeFactory);
        childFactory = factory;
        atomsPerGroup = atoms;
        configuration = config;
    }
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        childFactory.setSpecies(species);
    }
    
    /**
     * Constructs a new group.
     */
     public Atom makeAtom() {
        Atom group = newParentAtom();
        AtomTreeNodeGroup node = (AtomTreeNodeGroup)group.node;
        for(int i=0; i<atomsPerGroup; i++) {
            Atom childAtom = childFactory.makeAtom();
            childAtom.node.setParent(node);
        }
        configuration.initializePositions(node.childList);
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

     protected AtomFactory childFactory;
     private int atomsPerGroup;
     

}//end of AtomFactoryHomo
    
