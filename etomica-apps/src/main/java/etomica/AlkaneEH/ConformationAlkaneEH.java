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
	  protected final double CHbondL = 1.10/2;
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
		  
		  /////////////////////////// C0 and c[n-1] atoms ////////////////////////////////////////////////////
		  IAtom c0 = atomList.get(speciesAlkaneEH.indexC_3_1);// carbon atom in CH3 group on the left(label: "0")
		  c0.getPosition().E(new double[]{0.0,0.0,0.0});// put at the origin
		  
		  IAtom cEnd = atomList.get(speciesAlkaneEH.indexC_3_2);// carbon atom in CH3 group on the right(label:"n-1")
		  double cEndX = (numCarbons-1) * CCbondL * Math.sin(angleCCC_half);
		  double cEndY = (1-(numCarbons%2)) * CCbondL * Math.cos(angleCCC_half);// even number: not on x-axis; odd number: on x-axis
		  cEnd.getPosition().E(new double[]{cEndX,cEndY,0.0});
		 
		  ////////////////////// 6 H of CH3 groups  ////////////////////////////////////
		  IAtom h01 = atomList.get(numCarbons  * 3 );// on c0, in xOy plane
		  IAtom h02 = atomList.get(numCarbons );// on c0
		  IAtom h03 = atomList.get(numCarbons * 2);// on c0
		  IAtom hEnd1 = atomList.get(numCarbons * 3 + 1 );      // on c[n-1], in xOy plane
		  IAtom hEnd2 = atomList.get(numCarbons - 1 + numCarbons);// on c[n-1]
		  IAtom hEnd3 = atomList.get(numCarbons - 1 + numCarbons * 2 );//on c[n-1]
	  
		  double beta = angleCCH - angleCCC_half; // angle of (H-C-yAxis) 
		  double gamma = Math.PI- angleCCH;
		  double theta = Math.PI * 0.5 - angleCCC_half;
		  
		  ////////////////////// put 3H on the beginning C0 ////////////////////////////////////
		  
		  double c0M = CHbondL * Math.cos(gamma);// distance between C(H3) and center of HHH plane(M point)
		  Vector MPosition = space.makeVector();
		  Vector rNM = space.makeVector();
		  Vector NPosition = space.makeVector();
		  
		  double h01X = -CHbondL*Math.sin(beta);
		  double h01Y = CHbondL*Math.cos(beta);
		  h01.getPosition().E(new double[]{h01X,h01Y,0.0});

		  
		  MPosition.E(new double[]{ -c0M * Math.cos(theta), -c0M * Math.sin(theta), 0.0});
		  // get r(N) - r(M), direction: r(M)-r(H01), magnitude: |NM|=CHbond * sin*gamma / 2
		  rNM.Ev1Mv2(MPosition, h01.getPosition());
		  rNM.normalize();
		  rNM.TE(CHbondL * Math.sin(gamma)*0.5);
		  // r(N)=r(NM)+r(M)
		  NPosition.Ev1Pv2(rNM, MPosition);
		  
		  double h02X = NPosition.getX(0);
		  double h02Y = NPosition.getX(1);
		  double h02Z = CHbondL * Math.sin(gamma)*0.5 * Math.sqrt(3);
		  
		  h02.getPosition().E(new double[]{h02X, h02Y, h02Z});
		  h03.getPosition().E(new double[]{h02X, h02Y,-h02Z});
  
		  ////////////////////// put 3 H on the C[n-1] ////////////////////////////////////
		  double hEnd1X =  cEnd.getPosition().getX(0) + CHbondL * Math.sin(beta);
		  double hEnd1Y =  cEnd.getPosition().getX(1) - Math.pow(-1, numCarbons)* CHbondL * Math.cos(beta);
		  hEnd1.getPosition().E(new double[]{hEnd1X, hEnd1Y,0.0});
		    
		  double c1M1 = c0M;
		  Vector M1Position =  space.makeVector();
		  Vector rN1M1= space.makeVector();
		  Vector N1Position = space.makeVector();

		  double M1X = cEnd.getPosition().getX(0) + c1M1 * Math.cos(theta);
		  double M2X = cEnd.getPosition().getX(1) +Math.pow(-1, numCarbons)* c1M1 * Math.sin(theta);
		  M1Position.E(new double[]{M1X,M2X,0.0});
		  
		  // get r(N1) - r(M1), direction: r(M1)-r(H11), magnitude: |N1M1|=CHbond * sin*gamma / 2
		  rN1M1.Ev1Mv2(M1Position, hEnd1.getPosition());
		  rN1M1.normalize();
		  rN1M1.TE(CHbondL* Math.sin(gamma)*0.5);
		  // r(N1)=r(N1M1)+r(M1)
		  N1Position.Ev1Pv2(rN1M1, M1Position);
		  
		  double hEnd2X = N1Position.getX(0);
		  double hEnd2Y = N1Position.getX(1);
		  double hEnd2Z = CHbondL * Math.sin(gamma)*0.5 * Math.sqrt(3);
		  
		  hEnd2.getPosition().setX(0,hEnd2X);
		  hEnd2.getPosition().setX(1,hEnd2Y);
		  hEnd2.getPosition().setX(2,hEnd2Z);
	
		  hEnd3.getPosition().E(hEnd2.getPosition());
		  hEnd3.getPosition().setX(2, -hEnd2.getPosition().getX(2));
    	
		  /////////////////////////// C1~C[n-1] atoms ////////////////////////////////////////////////////
		  if (numCarbons > 2) {// true if the alkane chain has more than 2 carbons, need to put C(H2)
			  
			  for ( int m = 1; m < (numCH2+1); m++ ) {// put C of CH2 group in the chain, loop from 1st C(H2) to the last(numCH2)
				
				  IAtom c_ch2 = atomList.get(m);
				  double cX = m*CCbondL * Math.sin(angleCCC_half);
				  double cY = (m%2) * CCbondL * Math.cos(angleCCC_half);
				  
				  c_ch2.getPosition().E(new double[] { cX, cY, 0.0});

				  // put Ha and Hb on the chain, the first (numCarbons-2) number of H are specified first; then next (numCarbons-2) number of H put on the opposite direction
				  IAtom h_ch2_a = atomList.get( m + numCarbons );
				  IAtom h_ch2_b = atomList.get( m + numCarbons * 2 );

				  double h_ch2_aY = c_ch2.getPosition().getX(1) -Math.pow(-1, m)* CHbondL*Math.cos(angleHCH_half);
				  double h_ch2_aZ = CHbondL * Math.sin(angleHCH_half);
				  h_ch2_a.getPosition().E(new double[]{cX, h_ch2_aY,h_ch2_aZ});
				  ///////////////////// the other hydrogen//////////////////////////////////////
				  h_ch2_b.getPosition().E(h_ch2_a.getPosition());// the position is the same with n_ch2_a, except z-axis
				  h_ch2_b.getPosition().setX(2,-h_ch2_aZ);// change the sign of z-axis
				  
			  }	
		  }  
		  
		  
		  /*
		  System.out.println("................ Check initial conformation .................................");
          System.out.println("C0 position:" + c0.getPosition());
          System.out.println("C0's neighbor, C1 position:" + atomList.getAtom(1).getPosition());
          IVectorMutable vector_c0c1 = space.makeVector();
          vector_c0c1.Ev1Mv2(c0.getPosition(), atomList.getAtom(1).getPosition());
          double distance_c0c1 = Math.sqrt(vector_c0c1.squared());
          System.out.println("distance between C0 and C1, its neighbor, should be 1.535:" + distance_c0c1);
          
		  System.out.println("................ c0 and its H distance, should be 0.55");
          IVectorMutable vector_c0h01 = space.makeVector();
          vector_c0h01.Ev1Mv2(c0.getPosition(), atomList.getAtom(numCarbons).getPosition());
          double distance_c0h01 = Math.sqrt(vector_c0h01.squared());
          System.out.println("distance between C0 and h01, its neighbor:" + distance_c0h01);
          
          IVectorMutable vector_c0h02 = space.makeVector();
          vector_c0h02.Ev1Mv2(c0.getPosition(), atomList.getAtom(numCarbons*2).getPosition());
          double distance_c0h02 = Math.sqrt(vector_c0h02.squared());
          System.out.println("distance between C0 and h02, its neighbor:" + distance_c0h02);
          
          IVectorMutable vector_c0h03 = space.makeVector();
          vector_c0h03.Ev1Mv2(c0.getPosition(), atomList.getAtom(numCarbons*3).getPosition());
          double distance_c0h03 = Math.sqrt(vector_c0h03.squared());
          System.out.println("distance between C0 and h03, its neighbor:" + distance_c0h03);
          
          
          
		  System.out.println("................ HCH angle for CH3: ");
          IVectorMutable ch12 = space.makeVector();
          ch12.E(vector_c0h01);
          double cos12 = ch12.dot(vector_c0h02)/ Math.sqrt(vector_c0h01.squared())/Math.sqrt(vector_c0h02.squared());
          double angle12 = Math.acos(cos12) / Math.PI * 180.0;
          System.out.println("angle H1-C-H2 is: "+ angle12);
          
          IVectorMutable ch13 = space.makeVector();
          ch13.E(vector_c0h01);
          double cos13 = ch13.dot(vector_c0h03)/ Math.sqrt(vector_c0h01.squared())/Math.sqrt(vector_c0h03.squared());
          double angle13 = Math.acos(cos13) / Math.PI * 180.0;
          System.out.println("angle H1-C-H3 is: "+ angle13);
                             
          IVectorMutable ch23 = space.makeVector();
          ch23.E(vector_c0h02);
          double cos23 = ch23.dot(vector_c0h03)/ Math.sqrt(vector_c0h02.squared())/Math.sqrt(vector_c0h03.squared());
          double angle23 = Math.acos(cos23) / Math.PI * 180.0;
          System.out.println("angle H2-C-H3 is: "+ angle23);
          
          
          System.out.println("................ c1 and H01, H02, H03  distance, they should have the same value");
          IVectorMutable vector_c1h01 = space.makeVector();
          vector_c1h01.Ev1Mv2(atomList.getAtom(1).getPosition(), atomList.getAtom(numCarbons).getPosition());
          double distance_c1h01 = Math.sqrt(vector_c1h01.squared());
          System.out.println("distance between C1 and h01, its neighbor:" + distance_c1h01);
          
          IVectorMutable vector_c1h02 = space.makeVector();
          vector_c1h02.Ev1Mv2(atomList.getAtom(1).getPosition(), atomList.getAtom(numCarbons*2).getPosition());
          double distance_c1h02 = Math.sqrt(vector_c1h02.squared());
          System.out.println("distance between C1 and h02, its neighbor:" + distance_c1h02);
          
          IVectorMutable vector_c1h03 = space.makeVector();
          vector_c1h03.Ev1Mv2(atomList.getAtom(1).getPosition(), atomList.getAtom(numCarbons*3).getPosition());
          double distance_c1h03 = Math.sqrt(vector_c1h03.squared());
          System.out.println("distance between C1 and h03, its neighbor:" + distance_c1h03);
         
		       
          System.out.println("................ c2 and its H distance");
		  IVectorMutable vector_c2h1 = space.makeVector();
		  vector_c2h1.Ev1Mv2(atomList.getAtom(2).getPosition(), atomList.getAtom(numCarbons+2).getPosition());
          double distance_c2h1 = Math.sqrt(vector_c2h1.squared());
          System.out.println("distance between c2 and its h1:" + distance_c2h1);
          
          IVectorMutable vector_c2h2 = space.makeVector();
          vector_c2h2.Ev1Mv2(atomList.getAtom(2).getPosition(), atomList.getAtom(numCarbons*2+2).getPosition());
          double distance_c2h2 = Math.sqrt(vector_c2h2.squared());
          System.out.println("distance between C2 and its h2, its neighbor:" + distance_c2h2);
          
          
          System.out.println("................ HCH angle for C2H2, should be 107.8: ");
          IVectorMutable c2h1 = space.makeVector();
          c2h1.E(vector_c2h1);
          double cosHC2H = c2h1.dot(vector_c2h2)/ Math.sqrt(vector_c2h1.squared())/Math.sqrt(vector_c2h2.squared());
          double angleHC2H = Math.acos(cosHC2H) / Math.PI * 180.0;
          System.out.println("angle H1-C1-H2 is: "+ angleHC2H);
           
          
		  System.out.println("................ c[n-2] and c[n-1] distance");
		  IVectorMutable vector_ccEnd = space.makeVector();
		  vector_ccEnd.Ev1Mv2(atomList.getAtom(numCarbons-2).getPosition(), atomList.getAtom(numCarbons-1).getPosition());
          double distance_ccEnd = Math.sqrt(vector_ccEnd.squared());
          System.out.println("distance between C[n-2] and C[n-1] : :" + distance_ccEnd);
          
		  System.out.println("................ c[n-1] and its H distance");
          IVectorMutable vector_cEndh1 = space.makeVector();
          vector_cEndh1.Ev1Mv2(cEnd.getPosition(), atomList.getAtom(numCarbons*2-1).getPosition());
          double distance_cEndh1 = Math.sqrt(vector_cEndh1.squared());
          System.out.println("distance between C[n-1] and h1: " + distance_cEndh1);
          
          IVectorMutable vector_cEndh2 = space.makeVector();
          vector_cEndh2.Ev1Mv2(cEnd.getPosition(), atomList.getAtom(numCarbons*3-1).getPosition());
          double distance_c1h2 = Math.sqrt(vector_cEndh2.squared());
          System.out.println("distance between C[n-1] and h2 : " + distance_c1h2);
          
          IVectorMutable vector_c1h3 = space.makeVector();
          vector_c1h3.Ev1Mv2(cEnd.getPosition(), atomList.getAtom(numCarbons*3+1).getPosition());
          double distance_cEndh3 = Math.sqrt(vector_c1h3.squared());
          System.out.println("distance between C[n-1] and h3 : " + distance_cEndh3);
          
          
          System.out.println("................ c[n-2] and HEnd distance, should have the same value");
          IVectorMutable vector_cn2hEnd1 = space.makeVector();
          vector_cn2hEnd1.Ev1Mv2(atomList.getAtom(numCarbons-2).getPosition(), atomList.getAtom(numCarbons*2-1).getPosition());
          double distance_cn2hEnd1 = Math.sqrt(vector_cn2hEnd1.squared());
          System.out.println("distance between C[n-2] and h1: " + distance_cn2hEnd1);
          
          IVectorMutable vector_c2Endh2 = space.makeVector();
          vector_c2Endh2.Ev1Mv2(atomList.getAtom(numCarbons-2).getPosition(), atomList.getAtom(numCarbons*3-1).getPosition());
          double distance_cn2hEnd2 = Math.sqrt(vector_c2Endh2.squared());
          System.out.println("distance between C[n-2] and h2 : " + distance_cn2hEnd2);
          
          IVectorMutable vector_cn2hEnd3 = space.makeVector();
          vector_cn2hEnd3.Ev1Mv2(atomList.getAtom(numCarbons-2).getPosition(), atomList.getAtom(numCarbons*3+1).getPosition());
          double distance_cn2hEnd3 = Math.sqrt(vector_cn2hEnd3.squared());
          System.out.println("distance between C[n-2] and h3 : " + distance_cn2hEnd3);
          
          System.exit(0);
          
          */
	  }

}
