package etomica.atom;

import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.simulation.ISimulation;
import etomica.space.Space;
import etomica.species.Species;

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
    public AtomFactoryHomo(ISimulation sim, Species species) {
        this(sim, species, 1);
    }
    /**
     * @param space the coordinate factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param atoms the number of identical children per group (default is 1).
     */
    public AtomFactoryHomo(ISimulation sim, Species species, int atoms) {
        this(species, sim.getSpace(), atoms, new ConformationLinear(sim));
    }
    /**
     * @param space the coordinate factory
     * @param parentType the type instance of the atoms that are parents of those made by this factory
     * @param atoms the number of identical children per group (default is 1).
     * @param config the conformation applied to each group that is built (default is Linear).
     */
   public AtomFactoryHomo(Species species, Space space, int atoms, Conformation conformation) {
        super(new AtomTypeMolecule(species, new AtomPositionGeometricCenter(space)));
        atomsPerGroup = atoms;
        ((AtomTypeMolecule)atomType).setConformation(conformation);
    }
    
    /**
     * Constructs a new group.
     */
     public IAtom makeAtom() {
         isMutable = false;
         Molecule group = new Molecule(atomType);
         for(int i=0; i<atomsPerGroup; i++) {
             group.addChildAtom((IAtomLeaf)childFactory.makeAtom());
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
        this.childFactory = childFactory;
        ((AtomTypeMolecule)atomType).addChildType((AtomTypeLeaf)childFactory.getType());
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
    
