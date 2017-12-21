/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.AlkaneEH.SpeciesAlkaneEH;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.config.IConformation;
import etomica.space.Space;

/**
  *  Conformation for normal alkanes, explicit hydrogen
  *  Siepmann, TraPPE-EH
  * 
  * @author shu
  * 01-30-2013
  */
  public class ConformationAlkaneEH implements IConformation, java.io.Serializable{
	  protected final Space space;
	  private static final long serialVersionUID = 1L;

	  protected final double CCbondL = 1.535;
	  protected final double CHbondL = 1.10;
	  protected final double angleCCC =  112.7 * Math.PI/180; // CCC bond angle, in radians (mcWiggle will change this appropriately)
	  protected final double angleHCH =  107.8 * Math.PI/180; // HCH bond angle, in radians (mcWiggle will change this appropriately)
	  protected final double angleCCH = 110.7 * Math.PI/180;
	  protected double angleCCC_half = angleCCC / 2.0;
	  protected double angleHCH_half = angleHCH / 2.0;
	  protected SpeciesAlkaneEH speciesAlkaneEH ;
	  protected int numCarbons, numCH2, numH ; 
	  public ConformationAlkaneEH(Space space, SpeciesAlkaneEH speciesAlkaneEH){
		  this.space = space;
		  this.speciesAlkaneEH = speciesAlkaneEH;
		  this.numCarbons = speciesAlkaneEH.numCarbons;
		  this.numCH2 = speciesAlkaneEH.numCH2;
		  this.numH = speciesAlkaneEH.numH;
	  }

	  public void initializePositions(IAtomList atomList) {
		  
		  IAtom c0 = atomList.getAtom(speciesAlkaneEH.indexC_3_1);// carbon atom in CH3 group on the left
		  c0.getPosition().E(new double[] {0.0, 0.0,0.0});// put at the origin
		  IAtom c1 = atomList.getAtom(speciesAlkaneEH.indexC_3_2);// carbon atom in CH3 group on the right
		  c1.getPosition().E(new double[] { (numCarbons-1 ) * CCbondL * Math.sin(angleCCC_half), CCbondL * Math.cos(angleHCH_half), 0});
		  if (numCarbons %2 == 1) { // the right end is not on the x-axis
			  c1.getPosition().setX(1,0.0);//locate this C on the x-axis
		  }
		  // the last 6 H will be put near C(H3) on each end
		  IAtom h01 = atomList.getAtom(numCarbons  * 3 );// in xOy plane, special
		  IAtom h02 = atomList.getAtom(numCarbons );
		  IAtom h03 = atomList.getAtom(numCarbons * 2);
		  IAtom h11 = atomList.getAtom(numCarbons * 3 + 1 );// in xOy plane, special
		  IAtom h12 = atomList.getAtom(numCarbons - 1 + numCarbons);
		  IAtom h13 = atomList.getAtom(numCarbons - 1 + numCarbons * 2 );
	  
		  double angle = angleCCH - angleCCC_half; // angle of (H-C-yAxis)
		  h01.getPosition().E(new double[]{-CHbondL*Math.sin(angle), CHbondL*Math.cos(angle),0.0 });
		  h02.getPosition().E(new double[]{0.0, -CHbondL * Math.cos(angleHCH_half),CHbondL * Math.sin(angleHCH_half)});
		  h03.getPosition().E(new double[]{0.0, -CHbondL * Math.cos(angleHCH_half),-CHbondL * Math.sin(angleHCH_half)});
  
		  // assume c1 is not along the x-axis, thus the 2 Hs have larger y-axis values than the corresponding C
		  h11.getPosition().setX(0,c1.getPosition().getX(0) + CHbondL * Math.sin(angle));
		  h11.getPosition().setX(1, c1.getPosition().getX(1) - CHbondL * Math.cos(angle) );// this special H is smaller y-axis than C
		  h11.getPosition().setX(2,0.0);
    
		  h12.getPosition().setX(0,c1.getPosition().getX(0));
		  h12.getPosition().setX(1,c1.getPosition().getX(1)+CHbondL * Math.cos(angleHCH_half));
		  h12.getPosition().setX(2, CHbondL*Math.sin(angleHCH_half));
		  
		  if (numCarbons%2==1){
			  h11.getPosition().setX(1, c1.getPosition().getX(1)+CHbondL * Math.cos(angle));
			  h12.getPosition().setX(1,c1.getPosition().getX(1) - CHbondL * Math.cos(angleHCH_half));
		  }
		  h13.getPosition().E(h12.getPosition());
		  h13.getPosition().setX(2, -h12.getPosition().getX(2));
    	
		  
		  if (numCarbons > 2) {// true if the alkane chain has more than 2 carbons, nned to put C(H2)
			  for ( int m = 1; m < (numCH2+1); m++  ) {// put C of CH2 group in the chain, loop from 1st C(H2) to the last(numCH2)
				  IAtom c_ch2 = atomList.getAtom(m);
				  double c_xAxis = m*CCbondL * Math.sin(angleCCC_half);
				  c_ch2.getPosition().E(new double[] { c_xAxis, 0.0, 0.0});// assume m = even, C is on the x-axis

				  // put H on the chain, the first (numCarbons-2) number of H are specified first; then next (numCarbons-2) number of H put on the opposite direction
				  IAtom h_ch2_a = atomList.getAtom( m + numCarbons );
				  h_ch2_a.getPosition().setX(0,c_ch2.getPosition().getX(0));// x-axis = x-axis of Carbon
				  h_ch2_a.getPosition().setX(1,c_ch2.getPosition().getX(1) - CHbondL*Math.cos(angleHCH_half));
				  h_ch2_a.getPosition().setX(2,CHbondL * Math.sin(angleHCH_half)); // z-axis
				  if (m %2 == 1){// if C is not on the x-axis
					  double new_c_yAxis = CCbondL * Math.cos(angleCCC_half);
					  c_ch2.getPosition().setX(1, new_c_yAxis);// change c y-axis
					  h_ch2_a.getPosition().setX(1,new_c_yAxis + CHbondL*Math.cos(angleHCH_half));//change y-axis
				  }
        	
				  // the rest of the hydrogens
				  IAtom h_ch2_b = atomList.getAtom(m + numCarbons * 2 );
				  h_ch2_b.getPosition().E(h_ch2_a.getPosition());// the position is the same with n_ch2_a, except z-axis
				  h_ch2_b.getPosition().setX(2,-h_ch2_a.getPosition().getX(2));// change the sign of z-axis
       	
			  }	
	    
		  }       
	  }

}
