package etomica.models.nitrogen;

import etomica.api.IAtomTypeSphere;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomTypeSphere;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.space.ISpace;
import etomica.species.Species;

/**
 * 
 * 
 * Species nitrogen molecule (shell model) 
 * 
 * Reference: Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
 *             phases, JCP 112(15) 6745 (2000)
 *             
 * @author Tai Boon Tan
 *
 */
public class SpeciesN2ShellModel extends Species {

    public SpeciesN2ShellModel(ISpace space) {
        this(space, false);
    }
    
    public SpeciesN2ShellModel(ISpace space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;
        
        nType = new AtomTypeSphere(Nitrogen.INSTANCE, 3.1);
        pType = new AtomTypeSphere(new ElementSimple("P", 1.0), 0.0);
        addChildType(pType);
        addChildType(nType);
       
        setConformation(new ConformationNitrogenShellModel(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule nitrogenShellModel = new Molecule(this, 5);
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, nType) : new Atom(space, nType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, nType) : new Atom(space, nType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));


         conformation.initializePositions(nitrogenShellModel.getChildList());
         return nitrogenShellModel;
     }

     public IAtomTypeSphere getNitrogenType() {
         return nType;
     }

     public AtomTypeSphere getPType() {
         return pType;
     }


     public int getNumLeafAtoms() {
         return 5;
     }
    
     public final static int indexN1 = 0;
     public final static int indexN2 = 1;
     public final static int indexCenter = 2;
     public final static int indexP1left  = 3;
     public final static int indexP1right  = 4;


    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final boolean isDynamic;
    protected final AtomTypeSphere nType, pType;
}