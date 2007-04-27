package etomica.paracetamol;

import java.io.Serializable;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.normalmode.CoordinateDefinition;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class CoordinateDefinitionParacetamol extends CoordinateDefinition
        implements Serializable {

    public CoordinateDefinitionParacetamol(Space space, int numMolecule) {
        super(space.D() * numMolecule);
        this.space = space;
        work1 = space.makeVector();
        atomActionTranslateTo = new AtomActionTranslateTo(space);
        u = new double[coordinateDim];
    }

    /*
     * 
     */
    
    public double[] calcU(IAtom molecule) {
    	
    	AtomArrayList list = ((IAtomGroup)molecule).getChildList();
    	
    	int k =0;
    	
    	for (int i =0; i < list.size(); i++){
    		
    		IVector pos = ((AtomLeaf)((IAtomGroup)molecule).getChildList().get(i)).getPosition();
    		IVector site = getLatticePosition(((IAtomGroup)molecule).getChildList().get(i));
    		
    		for (int j =0; j <space.D(); j++){
    		
    			u[k] = pos.x(j) - site.x(j);
    			
    			k++;
    		}					
    	}		
        return u;
    
     }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(IAtom molecule) {
    	
    }

    public void setToU(IAtom molecule, double[] newU) {
    	
    	AtomArrayList list = ((IAtomGroup)molecule).getChildList();
    	
    	for (int i =0; i < list.size(); i++){
    		
    		IVector site = getLatticePosition(((IAtomGroup)molecule).getChildList().get(i));
    		
    		for (int j = 0; j < space.D(); j++) {
                work1.setX(j, site.x(j) + newU[j]);
            }
    		
    		atomActionTranslateTo.setDestination(work1);
            atomActionTranslateTo.actionPerformed(((IAtomGroup)molecule).getChildList().get(i));
    	}
     
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final IVector work1;
    protected final double[] u;
    protected final AtomActionTranslateTo atomActionTranslateTo;
}
