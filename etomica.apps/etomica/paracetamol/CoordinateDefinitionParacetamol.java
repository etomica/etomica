package etomica.paracetamol;

import java.io.Serializable;

import etomica.action.AtomActionTranslateTo;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtom;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space3d.Vector3D;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class CoordinateDefinitionParacetamol extends CoordinateDefinitionUpper
        implements Serializable {

    public CoordinateDefinitionParacetamol(Space space,int numCellMolecule, int atomLeaf) {
        super(space.D() * numCellMolecule * atomLeaf); // equals to coordinateDim
        this.space = space;
        work1 = space.makeVector();
        atomActionTranslateTo = new AtomActionTranslateTo(space);
        u = new double[coordinateDim];
        axes = new Vector3D [3];
        axes [0] = new Vector3D();
        axes [1] = new Vector3D();
        axes [2] = new Vector3D();
        com = new Vector3D();
        temp = new Vector3D();
    }

    /*
     * 
     */
    
    public double[] calcU(IAtom molecule) {
    	
    	AtomArrayList list = ((IAtomGroup)molecule).getChildList();
    	
    	int k =0;
    	
    	for (int i =0; i < list.size(); i++){
    		
    		IVector pos = ((IAtomPositioned)((IAtomGroup)molecule).getChildList().get(i)).getPosition();
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
    	
    	super.initNominalU(molecule);
    	
    	/*
    	 * Finding the vector for centre of mass
    	 */
    	AtomArrayList list = molecule.getParentGroup().getChildList();
    	for (int i= 0; i < list.size(); i++){
    		temp.Ea1Tv1(((AtomTypeLeaf)list.get(i).getType()).getMass(), ((IAtomPositioned)list.get(i)).getPosition()); //Need to input the MW of each molecule
    		com.PE(temp);
    	}
    	com.Ea1Tv1(1/151.2, com);
    	
    	
    	/*
    	 * Determine the Orientation of Each Molecule
    	 */
    	
    	IVector leafPos0 = ((IAtomPositioned)((IAtomGroup)molecule). 
				getParentGroup().getChildList().get(0)).getPosition();
    	IVector leafPos1 = ((IAtomPositioned)((IAtomGroup)molecule). 
				getParentGroup().getChildList().get(1)).getPosition();
    	IVector leafPos2 = ((IAtomPositioned)((IAtomGroup)molecule). 
				getParentGroup().getChildList().get(5)).getPosition();
    	
    	
    
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
    
    public void setCellManager(AtomAgentManager cell){
    	cellPosition = cell;
    }
    
    public void calcT(IVector k, double[] realT, double[] imaginaryT) {
        for (int i = 0; i < coordinateDim; i++) {
            realT[i] = 0;
            imaginaryT[i] = 0;
        }
        iterator.reset();
        // sum T over atoms
        
        
        int j = 0;
        for (IAtom atom = iterator.nextAtom(); atom != null;
             atom = iterator.nextAtom()) {
            double[] u0 = calcU(atom);
            
            
            IVector moleculeCellPosition = (IVector)cellPosition.getAgent(atom);
            
           
            double kR = k.dot(moleculeCellPosition);
            double coskR = Math.cos(kR);
            double sinkR = Math.sin(kR);
            for (int i = 0; i < u0.length; i++) {
                realT[j+i] += coskR * u0[i];
                imaginaryT[j+i] += sinkR * u0[i];
            }
            j+= u0.length;
            
            if(j == u0.length * 8){
            	j = 0;			//loop into the next cell of molecules
            }
        }

        for (int i = 0; i < coordinateDim; i++) {
            realT[i] /= sqrtN;
            imaginaryT[i] /= sqrtN;
        }

    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final IVector work1;
    protected final double[] u;
    protected final Vector3D [] axes;
    protected final Vector3D com, temp;
    protected final AtomActionTranslateTo atomActionTranslateTo;
    protected AtomAgentManager cellPosition;
}