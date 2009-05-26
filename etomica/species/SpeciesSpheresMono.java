package etomica.species;
import etomica.api.IAtom;
import etomica.api.IAtomType;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.AtomTypeSphere;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.space.ISpace;

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
    public SpeciesSpheresMono(ISimulation sim, ISpace _space) {
        this(sim, _space, new ElementSimple(sim));
    }
    
    public SpeciesSpheresMono(ISimulation sim, ISpace _space, Element element) {
        this(_space, sim.isDynamic(), new AtomTypeSphere(element));
    }
    
    public SpeciesSpheresMono(ISpace space, boolean isDynamic, AtomTypeSphere leafAtomType) {
        super(new AtomPositionFirstAtom());
        this.space = space;
        this.leafAtomType = leafAtomType;
        addChildType(leafAtomType);
        setConformation(new ConformationLinear(space, 1));
        this.isDynamic = isDynamic;
    }
    
    public IAtomType getLeafType() {
        return leafAtomType;
    }
    
    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         Molecule group = new Molecule(this, 1);
         group.addChildAtom(makeLeafAtom());
         return group;
     }

     protected IAtom makeLeafAtom() {
         return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                          : new Atom(space, leafAtomType);
     }

     public int getNumLeafAtoms() {
         return 1;
     }
     
     private static final long serialVersionUID = 1L;
     protected final ISpace space;
     protected final boolean isDynamic;
     protected final IAtomType leafAtomType;
}
