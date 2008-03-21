package etomica.species;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.AtomTypeSphere;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.space.Space;

/**
 * Species in which molecules are each made of a single spherical atom.
 * Does not permit multiatomic molecules.  The advantage of this species
 * over the multiatomic version (used with 1 atom), is that one layer of
 * the atom hierarchy is eliminated in SpeciesSpheresMono.  Each atom is
 * the direct child of the species agent (i.e., each atom is at the "molecule"
 * level in the hierarchy, without an intervening AtomGroup).
 * 
 * @author David Kofke
 */
public class SpeciesSpheresMono extends Species {

    /**
     * Constructs instance with a default element
     */
    public SpeciesSpheresMono(ISimulation sim, Space _space) {
        this(sim, _space, new ElementSimple(sim));
    }
    
    public SpeciesSpheresMono(ISimulation sim, Space _space, Element element) {
        this(_space, sim.isDynamic(), new AtomTypeSphere(element));
    }
    
    public SpeciesSpheresMono(Space space, boolean isDynamic, AtomTypeSphere leafAtomType) {
        super(new AtomPositionFirstAtom());
        this.space = space;
        this.leafAtomType = leafAtomType;
        addChildType(leafAtomType);
        setConformation(new ConformationLinear(space, 1));
        this.isDynamic = isDynamic;
    }
    
    public IAtomTypeLeaf getLeafType() {
        return leafAtomType;
    }
    
    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         Molecule group = new Molecule(this);
         group.addChildAtom(makeLeafAtom());
         return group;
     }

     protected IAtomLeaf makeLeafAtom() {
         return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                          : new AtomLeaf(space, leafAtomType);
     }

     public int getNumLeafAtoms() {
         return 1;
     }
     
     private static final long serialVersionUID = 1L;
     protected final Space space;
     protected final boolean isDynamic;
     protected final IAtomTypeLeaf leafAtomType;
}
