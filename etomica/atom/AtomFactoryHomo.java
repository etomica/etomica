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

/**
 * Builds an atom group that comprises a set of identically formed atoms or atom groups.
 * Construction of an instance must be followed by a call to setChildFactory, which identifies
 * the factory that makes the identically formed sub-atoms.
 * Default position definition is the geometric center (which is also the center of mass).
 *
 * @author David Kofke
 */

//child factory cannot be given in constructor because the AtomType the child factory will use cannot be
//constructed before this factory has made the AtomType used for its atom groups

public class AtomFactoryHomo extends AtomFactory {
    
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType) {
        this(space, sequencerFactory, parentType, 1);
    }
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param atoms the number of identical children per group (default is 1).
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType, int atoms) {
        this(space, sequencerFactory, parentType, atoms, new ConformationLinear(space));
    }
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param atoms the number of identical children per group (default is 1).
     * @param config the conformation applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
                            int atoms, Conformation config) {  
        this(space, sequencerFactory, parentType,
                  AtomTreeNodeGroup.FACTORY, atoms, config);
    }
 
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param nodeFactory makes nodes for each of the atoms built by this factory
     * @param atoms the number of identical children per group (default is 1).
     * @param config the conformation applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomSequencerFactory sequencerFactory, AtomTypeGroup parentType,
                            AtomTreeNodeFactory nodeFactory, int atoms, Conformation config) {
        super(space, new AtomTypeGroup(parentType, new AtomPositionGeometricCenter(space)), sequencerFactory, nodeFactory);
        atomsPerGroup = atoms;
        conformation = config;
    }
    
    public void setSpecies(Species species) {
        atomType.setSpecies(species);
        if (childFactory != null) {
            childFactory.setSpecies(species);
        }
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
        return group;
     }
     
    /**
     * Returns the subfactory that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory childFactory() {return childFactory;}

    /**
     * Sets the factory that makes the identical child atoms of an atom-group formed
     * by this factory.  This method should be called immediately after instantiation.
     * Subsequent attempts to invoke this method will throw an IllegalStateException.
     */
    public void setChildFactory(AtomFactory childFactory) {
        if (this.childFactory != null) throw new IllegalStateException("You can set the child factory only once!");
        if(childFactory == null) return;
        this.childFactory = childFactory;
        if(atomType.getSpecies() != null) childFactory.setSpecies(atomType.getSpecies());
    }
    
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
    
