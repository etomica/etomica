/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for TraPPE-Explicit hydrogen methane
 * C site and H site, NOTE the H is located in the middle of C-H bond 
 * @shu
 * 01-27-2013
 */
public class ConformationMethane implements IConformation {

    public ConformationMethane(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
    	double bond = 1.10 /2 ;
    	double angleEq =  109.5 * Math.PI/180;
    	double alpha = 30.0 * Math.PI/ 180;
    	double a = bond * Math.sin(Math.PI - angleEq ) * 2 * Math.cos(alpha) ;
    	// put H1 in the origin, H2 and H3 on the xy plane. C and H4 are located based on H1,H2 and H3
    	IAtom carbon = list.getAtom(SpeciesMethane.indexC);
    	IAtom h1 = list.getAtom(SpeciesMethane.indexH1);
    	IAtom h2 = list.getAtom(SpeciesMethane.indexH2);
    	IAtom h3 = list.getAtom(SpeciesMethane.indexH3);
    	IAtom h4 = list.getAtom(SpeciesMethane.indexH4);
    	double h1XAxis = -a * Math.sin(alpha) ; 
    	double h1YAxis = -a * Math.cos(alpha) ;
    	h1.getPosition().E(new double[] { 0.0 , 0.0 , 0.0 } ); 
    	h2.getPosition().E(new double[] {  h1XAxis, h1YAxis, 0.0 } ); 
    	h3.getPosition().E(new double[] { -h1XAxis, h1YAxis, 0.0 } ); 
    	
    	double carbonZAxis = bond * Math.cos( Math.PI - angleEq );
    	double carbonYAxis = - a / 2.0 / Math.cos(alpha) ; 
    	carbon.getPosition().E(new double[]{ 0.0,carbonYAxis, carbonZAxis }); 
    	h4.getPosition().E(new double[] { 0.0 , carbonYAxis, (carbonZAxis + bond) } ); 
    
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;

}
