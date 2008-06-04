package etomica.models.rowley;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IConformation;
import etomica.space.ISpace;

/**
 * Conformation for methanol as published in Rowley et al (2006)
 * 
 * OpenOffice Spreadsheet file, Methanol coordinates, used to calculate coordinates from bond lengths and angles: 
 * 
 * K.R. Schadel May 2008
 */
public class ConformationMethanol implements IConformation {

    public ConformationMethanol(ISpace space, boolean pointCharges) {
        this.space = space;
        this.pointCharges = pointCharges;
    }
    
    public void initializePositions(IAtomSet list){
        
    	// hydrogen attached to oxygen
        IAtomPositioned alpha_hydrogen = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexaH);
        alpha_hydrogen.getPosition().E(new double[] {0.0, 0.9114, 1.7172});
        
        IAtomPositioned oxygen = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexO);
        oxygen.getPosition().E(new double[] {0.0, 0.0, 1.4202});
                
        IAtomPositioned carbon = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexaC);
        carbon.getPosition().E(new double[] {0.0, 0.0, 0.0});
        
        // hydrogen on opposite of molecule from alpha hydrogen
        IAtomPositioned hydrogen_1 = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexH1);
        hydrogen_1.getPosition().E(new double[] {0.0, -1.0393, -0.3115});
        
        // hydrogen closer to alpha hydrogen
        IAtomPositioned hydrogen_2a = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexH2a);
        hydrogen_2a.getPosition().E(new double[] { 0.8920, 0.4850, -0.4101});
        
        // hydrogen closer to alpha hydrogen
        IAtomPositioned hydrogen_2b = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexH2b);
        hydrogen_2b.getPosition().E(new double[] {-0.8920, 0.4850, -0.4101});
        
        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        if(pointCharges) {
        	
        	// Methanol model with point charges at O, aC, aH
        	
        	// satellite site to model location of high electron density for oxygen
        	IAtomPositioned x = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexX);
        	x.getPosition().E(new double[] {0.0, -0.0753, 1.4749});	
        }
        else {

        	// Methanol model without point charges
        	
        	// satellite site to model location of high electron density for oxygen
        	IAtomPositioned x = (IAtomPositioned)list.getAtom(SpeciesMethanol.indexX);
        	x.getPosition().E(new double[] {0.0, -0.7774, 1.9845});
        }
       
        
        
    }
    
    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected boolean pointCharges;

}

