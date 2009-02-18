package etomica.models.traPPE;

import etomica.api.IAtomTypeSphere;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeSphere;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.ISpace;
import etomica.species.Species;

/**
 * Species for methanol with satellite site (Rowley et al 2006).
 */
public class SpeciesMethanol extends Species {

    public SpeciesMethanol(ISpace space) {
    	
        super(new AtomPositionGeometricCenter(space));
        
        this.space = space;
        
        cH3Type = new AtomTypeSphere(new ElementSimple("cH3", 1.0), 3.75); // diameter taken to be CH3-CH3 equilibrium LJ distance
        oType = new AtomTypeSphere(Oxygen.INSTANCE, 3.02); // diameter taken to be O-O equilibrium LJ distance
        hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0); // H-H equilibrium distance is not applicable 
        
        addChildType(cH3Type);
        addChildType(oType);
        addChildType(hType);
        
        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        setConformation(new ConformationMethanol(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule methanol = new Molecule(this);
         
         // The order in which the child atoms are added is important; it must match the site indices.
         methanol.addChildAtom(new Atom(space, cH3Type));
         methanol.addChildAtom(new Atom(space, oType));
         methanol.addChildAtom(new Atom(space, hType));
 
         conformation.initializePositions(methanol.getChildList());
         return methanol;
     }
     
     public IAtomTypeSphere getCH3Type() {
         return cH3Type;
     }
     
     public IAtomTypeSphere getOType() {
         return oType;
     }

     public IAtomTypeSphere getHType() {
         return hType;
     }

     public int getNumLeafAtoms() {
         return 3;
     }
    
    public final static int indexCH3  = 0;
    public final static int indexO   = 1;
    public final static int indexH  = 2;
  
    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final AtomTypeSphere cH3Type, oType, hType;
}
