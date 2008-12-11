package etomica.models.traPPE;

import etomica.api.IAtomList;
import etomica.api.IAtomPositioned;
import etomica.api.IConformation;
import etomica.space.ISpace;

/**
 * Conformation for TraPPE model of methanol: http://www.chem.umn.edu/groups/siepmann/trappe/molname.php#
 * 
 * CH3 group of methanol is united into a single site 
 * 
 * K.R. Schadel 2008
 */
public class ConformationMethanol implements IConformation {

    public ConformationMethanol(ISpace space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
    	
    	double bondCH3O = 1.43; // Angstroms
    	double bondOH = 0.95; // Angstroms
    	double angleEq = 108.50*Math.PI/180; // equilibrium bond angle in radians (mcWiggle will change this appropriately)
    	
    	IAtomPositioned cH3 = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexCH3);
        cH3.getPosition().E(new double[] {bondCH3O, 0.0, 0.0});
        
        IAtomPositioned oxygen = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexO);
        oxygen.getPosition().E(new double[] {0.0, 0.0, 0.0});
        
    	// hydrogen attached to oxygen
        IAtomPositioned hydrogen = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexH);
        hydrogen.getPosition().E(new double[] {bondOH*Math.cos(angleEq), bondOH*Math.sin(angleEq), 0.0});
   
        
    }
    
    private static final long serialVersionUID = 1L;
    protected final ISpace space;

}

