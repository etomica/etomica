package etomica.normalmode;

import etomica.api.IAtomType;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;
import etomica.species.Species;

/**
 * 
 * 
 * Species hard-sphere dimer molecule
 * 
 * @author Tai Boon Tan
 *
 */
public class SpeciesHSDimer extends Species {

    public SpeciesHSDimer(ISpace space) {
        this(space, false);
    }
    
    public SpeciesHSDimer(ISpace space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;
        
        dimerAtomType = new AtomTypeLeaf(new ElementSimple("P", 1.0));
        addChildType(dimerAtomType);

        setConformation(new ConformationHSDimer(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule hsDimer = new Molecule(this, 2);
         hsDimer.addChildAtom(isDynamic ? new AtomLeafDynamic(space, dimerAtomType) : new Atom(space, dimerAtomType));
         hsDimer.addChildAtom(isDynamic ? new AtomLeafDynamic(space, dimerAtomType) : new Atom(space, dimerAtomType));
         
         conformation.initializePositions(hsDimer.getChildList());
         return hsDimer;
     }

     public IAtomType getDimerAtomType() {
         return dimerAtomType;
     }

     public int getNumLeafAtoms() {
         return 2;
     }
    
    public final static int indexAtom1 = 0;
    public final static int indexAtom2 = 1;
    
    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final boolean isDynamic;
    protected final AtomTypeLeaf dimerAtomType;
}