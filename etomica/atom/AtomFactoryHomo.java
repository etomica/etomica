package etomica.atom;

import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * Builds an atom group that comprises a set of identically formed atoms or atom groups.
 * Construction of an instance must be followed by a call to setChildFactory, which identifies
 * the factory that makes the identically formed sub-atoms.
 * Default position definition is the geometric center (which is also the center of mass).
 *
 * @author David Kofke
 */
public class AtomFactoryHomo extends AtomFactory {
    
    /**
     * @param space the coordinate factory
     * @param sequencerFactory makes sequencers for each of the atoms built by this factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     */
    public AtomFactoryHomo(Simulation sim, AtomTypeGroup parentType) {
        this(sim, parentType, 1);
    }
    /**
     * @param space the coordinate factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param atoms the number of identical children per group (default is 1).
     */
    public AtomFactoryHomo(Simulation sim, AtomTypeGroup parentType, int atoms) {
        this(sim.space, parentType, atoms, new ConformationLinear(sim));
    }
    /**
     * @param space the coordinate factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param atoms the number of identical children per group (default is 1).
     * @param config the conformation applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomTypeGroup parentType,
                            int atoms, Conformation config) {  
        this(space, parentType, AtomTreeNodeGroup.FACTORY, atoms, config);
    }
 
    /**
     * @param space the coordinate factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param nodeFactory makes nodes for each of the atoms built by this factory
     * @param atoms the number of identical children per group (default is 1).
     * @param config the conformation applied to each group that is built (default is Linear).
     */
    public AtomFactoryHomo(Space space, AtomTypeGroup parentType,
                            AtomTreeNodeFactory nodeFactory, int atoms, Conformation config) {
        super(new AtomTypeGroup(new AtomPositionGeometricCenter(space)), nodeFactory);
        atomType.setParentType(parentType);
        atomsPerGroup = atoms;
        conformation = config;
    }
    
    /**
     * Constructs a new group.
     */
     public Atom makeAtom() {
         isMutable = false;
         Atom group = newParentAtom();
         AtomTreeNodeGroup node = (AtomTreeNodeGroup)group.getNode();
         for(int i=0; i<atomsPerGroup; i++) {
             Atom childAtom = childFactory.makeAtom();
             childAtom.getNode().setParent(node);
         }
         return group;
     }
     
    /**
     * Returns the subfactory that produces each of the identical atoms
     * in the group made by this factory.
     */
    public AtomFactory getChildFactory() {
        return childFactory;
    }

    /**
     * Sets the factory that makes the identical child atoms of an atom-group formed
     * by this factory.  This method should be called immediately after instantiation.
     * Subsequent attempts to invoke this method will throw an IllegalStateException.
     * 
     * @throws IllegalStateException if invoked more than once for an instance
     */
    public void setChildFactory(AtomFactory childFactory) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        if (this.childFactory != null) throw new IllegalStateException("You can set the child factory only once!");
        if(childFactory == null) return;
        this.childFactory = childFactory;
    }
    
    /**
     * Specifies the number of child atoms in each atom constructed by this factory.
     * 
     * @param na The new number of atoms per group
     */
    public void setNumChildAtoms(int na) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        atomsPerGroup = na;
    }

     public int getNumChildAtoms() {
         return atomsPerGroup;
     }

     public int getNumTreeAtoms() {
         return 1 + atomsPerGroup*childFactory.getNumTreeAtoms();
     }
     
     public int getNumLeafAtoms() {
         return atomsPerGroup*childFactory.getNumLeafAtoms();
     }
     
     private static final long serialVersionUID = 1L;
     protected AtomFactory childFactory;
     protected int atomsPerGroup;
}
    
