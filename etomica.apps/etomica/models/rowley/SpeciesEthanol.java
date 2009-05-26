package etomica.models.rowley;


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
 * Species for ethanol with satellite site (Rowley et al 2006).
 */
public class SpeciesEthanol extends Species {

    public SpeciesEthanol(ISpace space, boolean pointCharges) {
    	
        super();
        
        this.space = space;
        
        oType = new AtomTypeSphere(Oxygen.INSTANCE, 2.0); // diameter NOT taken to be O-O equilibrium distance
        acType = new AtomTypeSphere(Carbon.INSTANCE, 2.0); // diameter NOT taken to be aC-aC equilibrium distance
        ahType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0); // diameter NOT taken to be aH-aH equilibrium distance
        cType = new AtomTypeSphere(Carbon.INSTANCE, 2.0); // diameter NOT taken to be aC-aC equilibrium distance
        hType = new AtomTypeSphere(Hydrogen.INSTANCE, 2.0); // diameter NOT taken to be H-H equilibrium distance
        xType = new AtomTypeSphere(new ElementSimple("X", 1.0), 2.0); // diameter NOT taken to be X-X equilibrium distance
        
        addChildType(oType);
        addChildType(acType);
        addChildType(ahType);
        addChildType(cType);
        addChildType(hType);
        addChildType(xType);
        
        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        setConformation(new ConformationEthanol(space, pointCharges)); 
     }

     public IMolecule makeMolecule() {
         Molecule ethanol = new Molecule(this, 10);
         
         // The order in which the child atoms are added is important; it must match the order of site indices below.
         ethanol.addChildAtom(new Atom(space, oType));
         ethanol.addChildAtom(new Atom(space, acType));
         ethanol.addChildAtom(new Atom(space, ahType));
         ethanol.addChildAtom(new Atom(space, cType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, hType));
         ethanol.addChildAtom(new Atom(space, xType));
         conformation.initializePositions(ethanol.getChildList());
         return ethanol;
     }
     
     public IAtomTypeSphere getOxygenType() {
         return oType;
     }
     
     public IAtomTypeSphere getAlphaCarbonType() {
         return acType;
     }
     
     public IAtomTypeSphere getAlphaHydrogenType() {
         return ahType;
     }
     
     public IAtomTypeSphere getCarbonType() {
         return cType;
     }

     public IAtomTypeSphere getHydrogenType() {
         return hType;
     }

     public IAtomTypeSphere getXType() {
         return xType;
     }

     public int getNumLeafAtoms() {
         return 10;
     }
    
    public final static int indexO   = 0;
    public final static int indexaC  = 1;
    public final static int indexaH  = 2; // ahType
    public final static int indexC   = 3;
    public final static int indexH1a = 4; // hType
    public final static int indexH1b = 5; // hType
    public final static int indexH2a = 6; // hType
    public final static int indexH2b = 7; // hType
    public final static int indexH2c = 8; // hType
    public final static int indexX   = 9;

    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final AtomTypeSphere oType, acType, ahType, cType, hType, xType;
}
