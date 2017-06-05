/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.AlkaneEH;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.config.IConformation;
import etomica.space.Space;

/**
 * Conformation for TraPPE-Explicit hydrogen methane
 * C site and H site, NOTE the H is located in the middle of C-H bond 
 * @author shu
 * 01-27-2013
 */
public class ConformationMethane implements IConformation {

    public ConformationMethane(Space space) {
        this.space = space;
    }
    
    public void initializePositions(IAtomList list){
    	double bond = 1.10 /2;
   // 	double angleEq =  109.5 * Math.PI/180;
    	
    	double angleEq =  2 * Math.acos(1.0/Math.sqrt(3.0));
    	//System.out.println("angleEq:"+angleEq);
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
    	if (false){    	  
	    	System.out.println("................ Check initial conformation .................................");
	    	System.out.println("C position:" + carbon.getPosition());
	             
	    	System.out.println("................ c and its H1 distance:");
	        Vector vector_ch1 = space.makeVector();
	        vector_ch1.Ev1Mv2(carbon.getPosition(), h1.getPosition());
	        double distance_ch1 = Math.sqrt(vector_ch1.squared());
	    	System.out.println("h1 position:" + h1.getPosition());
	        System.out.println("distance between C and h1, its neighbor:" + distance_ch1);
	
	    	System.out.println("................ c and its H2 distance:");
	        Vector vector_ch2 = space.makeVector();
	        vector_ch2.Ev1Mv2(carbon.getPosition(), h2.getPosition());
	        double distance_ch2 = Math.sqrt(vector_ch2.squared());
	    	System.out.println("h2 position:" + h2.getPosition());
	        System.out.println("distance between C and h2, its neighbor:" + distance_ch2);
	       
	    	System.out.println("................ c and its H3 distance:");
	        Vector vector_ch3 = space.makeVector();
	        vector_ch3.Ev1Mv2(carbon.getPosition(), h3.getPosition());
	        double distance_ch3 = Math.sqrt(vector_ch3.squared());
	    	System.out.println("h3 position:" + h3.getPosition());
	        System.out.println("distance between C and h3, its neighbor:" + distance_ch3);
	        
	    	System.out.println("................ c and its H4 distance:");
	        Vector vector_ch4 = space.makeVector();
	        vector_ch4.Ev1Mv2(carbon.getPosition(), h4.getPosition());
	        double distance_ch4 = Math.sqrt(vector_ch4.squared());
	    	System.out.println("h4 position:" + h4.getPosition());
	        System.out.println("distance between C and h4, its neighbor:" + distance_ch4);
	        
	        System.out.println("................ HCH angle : ");
	        Vector ch12 = space.makeVector();
	        ch12.E(vector_ch1);
	        double cos12 = ch12.dot(vector_ch2)/ Math.sqrt(vector_ch1.squared())/Math.sqrt(vector_ch2.squared());
	        double angle1C2 = Math.acos(cos12) / Math.PI * 180.0;
	        System.out.println("angle H1-C-H2 is: "+ angle1C2);
	        
	        Vector ch13 = space.makeVector();
	        ch13.E(vector_ch1);
	        double cos13 = ch13.dot(vector_ch3)/ Math.sqrt(vector_ch1.squared())/Math.sqrt(vector_ch3.squared());
	        double angle1C3 = Math.acos(cos13) / Math.PI * 180.0;
	        System.out.println("angle H1-C-H3 is: "+ angle1C3);
	        
	        Vector ch14 = space.makeVector();
	        ch14.E(vector_ch1);
	        double cos14 = ch14.dot(vector_ch4)/ Math.sqrt(vector_ch1.squared())/Math.sqrt(vector_ch4.squared());
	        double angle1C4 = Math.acos(cos14) / Math.PI * 180.0;
	        System.out.println("angle H1-C-H4 is: "+ angle1C4);
	                                            
	        Vector ch23 = space.makeVector();
	        ch23.E(vector_ch2);
	        double cos23 = ch23.dot(vector_ch3)/ Math.sqrt(vector_ch2.squared())/Math.sqrt(vector_ch3.squared());
	        double angle2C3 = Math.acos(cos23) / Math.PI * 180.0;
	        System.out.println("angle H2-C-H3 is: "+ angle2C3);
	        
	        Vector ch24 = space.makeVector();
	        ch24.E(vector_ch2);
	        double cos24 = ch24.dot(vector_ch4)/ Math.sqrt(vector_ch2.squared())/Math.sqrt(vector_ch4.squared());
	        double angle2C4 = Math.acos(cos24) / Math.PI * 180.0;
	        System.out.println("angle H2-C-H4 is: "+ angle2C4);
	        
	        Vector ch34 = space.makeVector();
	        ch34.E(vector_ch3);
	        double cos34 = ch34.dot(vector_ch4)/ Math.sqrt(vector_ch3.squared())/Math.sqrt(vector_ch4.squared());
	        double angle3C4 = Math.acos(cos34) / Math.PI * 180.0;
	        System.out.println("angle H3-C-H4 is: "+ angle3C4);
	        System.exit(0);
    	}
    }
    
    private static final long serialVersionUID = 1L;
    protected final Space space;

}
