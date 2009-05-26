package etomica.species;
import etomica.api.IAtom;
import etomica.api.IAtomType;
import etomica.api.IConformation;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeSphere;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.space.ISpace;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */
public class SpeciesSpheres extends Species {

    public SpeciesSpheres(ISimulation sim, ISpace _space) {
        this(sim, _space, 1);
    }
    public SpeciesSpheres(ISimulation sim, ISpace _space, int nA) {
        this(sim, _space, nA, new ElementSimple(sim));
    }
    
    public SpeciesSpheres(ISimulation sim, ISpace _space, int nA, Element leafElement) {
        this(sim, nA, leafElement, new ConformationLinear(_space), _space);
    }
    
    public SpeciesSpheres(ISimulation sim, int nA, Element leafElement,
    		              IConformation conformation, ISpace _space) {
        this(_space, sim.isDynamic(), nA, new AtomTypeSphere(leafElement), conformation);
    }
    
    public SpeciesSpheres(ISpace _space, boolean isDynamic, int nA, IAtomType leafAtomType, IConformation conformation) {
        super(new AtomPositionGeometricCenter(_space));
        this.space = _space;
        addChildType(leafAtomType);
        setNumLeafAtoms(nA);
        setConformation(conformation);
        this.leafAtomType = leafAtomType;
        this.isDynamic = isDynamic;
    }
    
    public IAtomType getLeafType() {
        return getAtomType(0);
    }
    
    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         Molecule group = new Molecule(this, atomsPerGroup);
         for(int i=0; i<atomsPerGroup; i++) {
             group.addChildAtom(makeLeafAtom());
         }
         conformation.initializePositions(group.getChildList());
         return group;
     }
     
     protected IAtom makeLeafAtom() {
         return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                          : new Atom(space, leafAtomType);
     }

    /**
     * Specifies the number of child atoms in each atom constructed by this factory.
     * 
     * @param na The new number of atoms per group
     */
    public void setNumLeafAtoms(int na) {
        atomsPerGroup = na;
    }

     public int getNumLeafAtoms() {
         return atomsPerGroup;
     }

     private static final long serialVersionUID = 1L;
     protected final boolean isDynamic;
     protected final ISpace space;
     protected int atomsPerGroup;
     protected final IAtomType leafAtomType;
}
