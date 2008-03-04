package etomica.species;
import java.lang.reflect.Constructor;

import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.Conformation;
import etomica.config.ConformationLinear;
import etomica.space.Space;

/**
 * Species in which molecules are made of arbitrary number of spheres,
 * with each sphere having the same mass and size (same type).
 * 
 * @author David Kofke
 */

public class SpeciesSpheres extends Species {

    public SpeciesSpheres(ISimulation sim) {
        this(sim, 1);
    }
    public SpeciesSpheres(ISimulation sim, int nA) {
        this(sim, nA, new ElementSimple(sim));
    }
    
    public SpeciesSpheres(ISimulation sim, int nA, Element leafElement) {
        this(sim, nA, leafElement, new ConformationLinear(sim));
    }
    
    public SpeciesSpheres(ISimulation sim, int nA, Element leafElement, Conformation conformation) {
        this(sim.getSpace(), sim.isDynamic(), nA, new AtomTypeSphere(leafElement), conformation);
    }
    
    public SpeciesSpheres(Space space, boolean isDynamic, int nA, AtomTypeLeaf leafAtomType, Conformation conformation) {
        super(new AtomTypeMolecule(new AtomPositionGeometricCenter(space)));
        this.space = space;
        atomType.addChildType(leafAtomType);
        setNumLeafAtoms(nA);
        atomType.setConformation(conformation);
        this.leafAtomType = leafAtomType;
        this.isDynamic = isDynamic;
    }
    
    public AtomTypeLeaf getLeafType() {
        return getMoleculeType().getChildTypes()[0];
    }
    
    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         isMutable = false;
         Molecule group = new Molecule(atomType);
         for(int i=0; i<atomsPerGroup; i++) {
             group.addChildAtom(makeLeafAtom());
         }
         return group;
     }
     
     protected IAtomLeaf makeLeafAtom() {
         return isDynamic ? new AtomLeafDynamic(space, leafAtomType)
                          : new AtomLeaf(space, leafAtomType);
     }

    /**
     * Specifies the number of child atoms in each atom constructed by this factory.
     * 
     * @param na The new number of atoms per group
     */
    public void setNumLeafAtoms(int na) {
        if (!isMutable) {
            throw new IllegalStateException("Factory is not mutable");
        }
        atomsPerGroup = na;
    }

     public int getNumLeafAtoms() {
         return atomsPerGroup;
     }

     public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class, Integer.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(constructor,new Object[]{new Integer(atomsPerGroup)});
    }
    
     private static final long serialVersionUID = 1L;
     protected final boolean isDynamic;
     protected final Space space;
     protected int atomsPerGroup;
     protected final AtomTypeLeaf leafAtomType;
}
