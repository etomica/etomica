package etomica.paracetamol;

import java.io.Serializable;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.CoordinateDefinitionMolecule;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * CoordinateDefinition implementation for molecules. The class takes the first
 * space.D values of u to be real space displacements of the molecule center of
 * mass from its nominal position. Subclasses should add additional u values for
 * intramolecular degrees of freedom.
 * 
 * @author Andrew Schultz & Tai Tan
 */
public class CoordinateDefinitionParacetamol extends CoordinateDefinitionMolecule
        implements Serializable {

    public CoordinateDefinitionParacetamol(Phase phase, Primitive primitive, Basis basis) {
        super(phase, primitive, 3, basis);
        axes = new IVector [3];
        axes [0] = phase.getSpace().makeVector();
        axes [1] = phase.getSpace().makeVector();
        axes [2] = phase.getSpace().makeVector();
        com = phase.getSpace().makeVector();
        temp = phase.getSpace().makeVector();
    }

    /*
     * 
     */
    
    public double[] calcU(AtomSet molecules) {
        super.calcU(molecules);
        
        //XXX need to calculate u components related to orientation
        
        return u;
     }

    /**
     * Override if nominal U is more than the lattice position of the molecule
     */
    public void initNominalU(AtomSet molecules) {
    	//XXX need to save info about orientation
        
        //XXX loop over molecules
        
        IAtomGroup molecule = null;
        
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

    public void setToU(AtomSet molecules, double[] newU) {
        super.setToU(molecules, newU);
        
        //XXX need to set orientation of molecules
    }

    private static final long serialVersionUID = 2L;
    protected final IVector [] axes;
    protected final IVector com, temp;
}
