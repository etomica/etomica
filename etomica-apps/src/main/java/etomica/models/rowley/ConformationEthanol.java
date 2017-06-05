/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.atom.IAtomList;
import etomica.atom.IAtom;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for ethanol (positioned in standard orientation) as published in Rowley et al (2006)
 * 
 * K.R. Schadel 2008
 */
public class ConformationEthanol implements IConformation {

    public ConformationEthanol(Space space, boolean pointCharges) {
        this.space = space;
        this.pointCharges = pointCharges;
    }
    
    public void initializePositions(IAtomList list){
    	
    	/*
    	 *  Standard orientation for ethanol molecule sites: 
    	 *  	The alpha carbon (carbon bonded to oxygen) is at the origin.
    	 *  	The bond between oxygen and alpha carbon is on the positive z-axis
    	 *  	The bond between oxygen and alpha hydrogen is in the xz plane (positive x-axis part).
    	 * 
    	 *  The standard coordinates in eth01.out provided by R. Rowley do NOT correspond to the standard orientation described in his 2006 paper.
    	 *  The calculation of the standard-orientation coordinates from these coordinates are provided in Coordinates.ods.
    	 */
 
 
    	
    	IAtom oxygen = list.getAtom(SpeciesEthanol.indexO);
        oxygen.getPosition().E(new double[]         {0.000000,    0.000000,   1.425003});
        
        // carbon attached to oxygen
        IAtom alpha_carbon = list.getAtom(SpeciesEthanol.indexaC);
        alpha_carbon.getPosition().E(new double[]   { 0.000000,    0.000000,   0.000000});
        
        // carbon NOT attached to oxygen
        IAtom carbon = list.getAtom(SpeciesEthanol.indexC);
        carbon.getPosition().E(new double[]         {-1.439886,    0.000000,  -0.449889});
        
        // hydrogen attached to oxygen
        IAtom alpha_hydrogen = list.getAtom(SpeciesEthanol.indexaH);
        alpha_hydrogen.getPosition().E(new double[] { 0.910604,    0.000000,   1.727510});
        
        // hydrogen attached to alpha carbon
        IAtom hydrogen_1a = list.getAtom(SpeciesEthanol.indexH1a);
        hydrogen_1a.getPosition().E(new double[]    { 0.515370,   -0.883917,   -0.383467});
        
        // hydrogen attached to alpha carbon
        IAtom hydrogen_1b = list.getAtom(SpeciesEthanol.indexH1b);
        hydrogen_1b.getPosition().E(new double[]    { 0.515370,    0.883917,   -0.383467});
        
        // hydrogen attached to non-alpha carbon
        IAtom hydrogen_2a = list.getAtom(SpeciesEthanol.indexH2a);
        hydrogen_2a.getPosition().E(new double[]    {-1.952009,    0.882056,    -0.072720});
        
        // hydrogen attached to non-alpha carbon
        IAtom hydrogen_2b = list.getAtom(SpeciesEthanol.indexH2b);
        hydrogen_2b.getPosition().E(new double[]    {-1.952009,   -0.882056,    -0.072720});
        
        // hydrogen attached to non-alpha carbon
        IAtom hydrogen_2c = list.getAtom(SpeciesEthanol.indexH2c);
        hydrogen_2c.getPosition().E(new double[]    {-1.500395,    0.000000,    -1.536664});
        
        // The satellite site, X, is closer to the oxygen atom in the model with point charges.
        // Calculation of X position made in Coordinates.ods
        if(pointCharges) {
        	
        	// Ethanol model with point charges at O, aC, aH
        	
        	// satellite site to model location of high electron density for oxygen
        	IAtom x = list.getAtom(SpeciesEthanol.indexX);
        	x.getPosition().E(new double[]          {-0.274574,    0.000000,    1.623110});
        }
        else {

        	// Ethanol model without point charges
        	
        	// satellite site to model location of high electron density for oxygen
        	IAtom x = list.getAtom(SpeciesEthanol.indexX);
        	x.getPosition().E(new double[]          {-0.898492,    0.000000,    2.073284});
        }
       
        
        
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected boolean pointCharges;

}


/*

Coordinates from etho1.out

1          8             0        0.000000    2.000100   -1.041997  O
2          1             0        0.000000    2.302608   -1.952601	aH
3          6             0        0.000000    3.141836   -0.189312	aC
4          6             0        0.000000    2.640704    1.233551	C
5          1             0       -0.883917    3.757460   -0.372778 	H (H1a)
6          1             0        0.883917    3.757460   -0.372778  H (H1b)
7          1             0        0.882056    2.032069    1.418184  H (H2a)
8          1             0       -0.882056    2.032069    1.418184  H (H2b)
9          1             0        0.000000    3.475239    1.932330  H (H2c)
*/
