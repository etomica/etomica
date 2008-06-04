package etomica.models.rowley;


import etomica.api.IAtomTypeSphere;
import etomica.api.IMolecule;
import etomica.atom.AtomLeaf;
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

    public SpeciesMethanol(ISpace space, boolean pointCharges) {
    	
        super(new AtomPositionGeometricCenter(space));
        
        this.space = space;
        
        oType = new AtomTypeSphere(Oxygen.INSTANCE, 2.0); // diameter NOT taken to be O-O equilibrium distance
        cType = new AtomTypeSphere(Carbon.INSTANCE, 2.0); // diameter NOT taken to be aC-aC equilibrium distance
        ahType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0); // diameter NOT taken to be aH-aH equilibrium distance
        hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0); // diameter NOT taken to be H-H equilibrium distance
        xType = new AtomTypeSphere(new ElementSimple("X", 1.0), 2.0); // diameter NOT taken to be X-X equilibrium distance
        
        addChildType(oType);
        addChildType(cType);
        addChildType(ahType);
        addChildType(hType);
        addChildType(xType);
        
        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        setConformation(new ConformationMethanol(space, pointCharges)); 
     }

     public IMolecule makeMolecule() {
         Molecule methanol = new Molecule(this);
         
         // The order in which the child atoms are added is important; it must match the site indices.
         methanol.addChildAtom(new AtomLeaf(space, oType));
         methanol.addChildAtom(new AtomLeaf(space, cType));
         methanol.addChildAtom(new AtomLeaf(space, ahType));
         methanol.addChildAtom(new AtomLeaf(space, hType));
         methanol.addChildAtom(new AtomLeaf(space, hType));
         methanol.addChildAtom(new AtomLeaf(space, hType));
         methanol.addChildAtom(new AtomLeaf(space, xType));
         conformation.initializePositions(methanol.getChildList());
         return methanol;
     }
     
     public IAtomTypeSphere getOxygenType() {
         return oType;
     }
     
     public IAtomTypeSphere getCarbonType() {
         return cType;
     }

     public IAtomTypeSphere getAlphaHydrogenType() {
         return ahType;
     }

     public IAtomTypeSphere getHydrogenType() {
         return hType;
     }

     public IAtomTypeSphere getXType() {
         return xType;
     }

     public int getNumLeafAtoms() {
         return 7;
     }
    
    public final static int indexO   = 0;
    public final static int indexaC  = 1;
    public final static int indexaH  = 2; // ahType
    public final static int indexH1  = 3; // hType
    public final static int indexH2a = 4; // hType
    public final static int indexH2b = 5; // hType
    public final static int indexX   = 6;

    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final AtomTypeSphere oType, cType, ahType, hType, xType;
}
