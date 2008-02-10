package etomica.species;
import java.lang.reflect.Constructor;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionFirstAtom;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomLeaf;
import etomica.atom.IMolecule;
import etomica.atom.Molecule;
import etomica.chem.elements.Element;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConformationLinear;
import etomica.simulation.ISimulation;
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
    public SpeciesSpheresMono(ISimulation sim) {
        this(sim, new ElementSimple(sim));
    }
    
    public SpeciesSpheresMono(ISimulation sim, Element element) {
        this(sim.getSpace(), sim.isDynamic(), new AtomTypeSphere(element));
    }
    
    public SpeciesSpheresMono(Space space, boolean isDynamic, AtomTypeSphere leafAtomType) {
        super(new AtomTypeMolecule(new AtomPositionFirstAtom()));
        this.space = space;
        this.leafAtomType = leafAtomType;
        atomType.addChildType(leafAtomType);
        atomType.setConformation(new ConformationLinear(space, 1));
        this.isDynamic = isDynamic;
    }
    
    public AtomTypeLeaf getLeafType() {
        return leafAtomType;
    }
    
    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         isMutable = false;
         Molecule group = new Molecule(atomType);
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
     
     public SpeciesSignature getSpeciesSignature() {
        Constructor constructor = null;
        try {
            constructor = this.getClass().getConstructor(new Class[]{ISimulation.class,Element.class});
        }
        catch(NoSuchMethodException e) {
            System.err.println("you have no constructor.  be afraid");
        }
        return new SpeciesSignature(constructor,new Object[]{getLeafType().getElement()});
    }
    
     private static final long serialVersionUID = 1L;
     protected final Space space;
     protected final boolean isDynamic;
     protected final AtomTypeLeaf leafAtomType;
}
