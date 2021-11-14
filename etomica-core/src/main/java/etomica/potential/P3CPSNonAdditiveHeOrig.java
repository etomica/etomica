/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.dimensions.Angle;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;

/**
 * Ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009.
 * 
 * This class is very slow.  Please do not use it.
 * @author kate
 */
public class P3CPSNonAdditiveHeOrig implements PotentialSoft {

    public P3CPSNonAdditiveHeOrig(Space space) {
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
        setAngle(Math.PI);
        gradient = new Vector[3];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        gradient[2] = space.makeVector();
        
        
    }

    public double energy(IAtomList atomSet) {
    	
    	setA();
        setAlpha();
        setBeta3();
        setZ3();
        setBeta4220();
        setBeta4211();
        setZ4220();
        setZ4211();
        
        //Operate on duplicate of atomSet
        
        IAtomList atomSet2 = atomSet;
        
        IAtom atomA = atomSet2.get(0);
        IAtom atomB = atomSet2.get(1);
        IAtom atomC = atomSet2.get(2);
        
       
        drAB.Ev1Mv2(atomA.getPosition(),atomB.getPosition());
        drAC.Ev1Mv2(atomA.getPosition(),atomC.getPosition());
        drBC.Ev1Mv2(atomB.getPosition(),atomC.getPosition());
        
        drAB.TE(1.0/AngstromPerBohrRadius);
        drAC.TE(1.0/AngstromPerBohrRadius);
        drBC.TE(1.0/AngstromPerBohrRadius);
        
        
        double RAB = Math.sqrt(drAB.squared());
        double RAC = Math.sqrt(drAC.squared());
        double RBC = Math.sqrt(drBC.squared());
        
        //Supplementary 0 energy region 
        
        if (RAB<2.23 || RAC<2.23 || RBC<2.23) {
        	return 0.0;
        }
        
        double costhetaA =  drAB.dot(drAC)/(RAB*RAC);
        double costhetaB = -drAB.dot(drBC)/(RAB*RBC);
        double costhetaC =  drAC.dot(drBC)/(RAC*RBC);
        
        
        
        double thetaA; double thetaB; double thetaC;
        
        if (costhetaA > 1) { thetaA = 0;}
        else if (costhetaA < -1) { thetaA = Math.PI;}
        else {thetaA = Math.acos(costhetaA);}
        
        if (costhetaB > 1) {thetaB = 0;}
        else if (costhetaB < -1) {thetaB = Math.PI;}
        else {thetaB = Math.acos(costhetaB);}
        
        if (costhetaC > 1) {thetaC = 0;}
        else if (costhetaC < -1) {thetaC = Math.PI;}
        else {thetaC = Math.acos(costhetaC);}
        
        //System.out.println(thetaA + " " + thetaB + " " + thetaC + "  " + (thetaA+thetaB+thetaC)/Math.PI);
        
        double Vexp = 0;

        for (int k3 = 0; k3<=4; k3++) {
        	for (int k2 = 0; k2<=k3; k2++) {
        		for (int k1 = 0; k1<=k2; k1++) {
        			
        			double P = legendreP(k1,costhetaA)*legendreP(k2,costhetaB)*legendreP(k3,costhetaC);
        			P =    P + legendreP(k1,costhetaA)*legendreP(k2,costhetaC)*legendreP(k3,costhetaB);
        			P =    P + legendreP(k1,costhetaB)*legendreP(k2,costhetaA)*legendreP(k3,costhetaC);
        			P =    P + legendreP(k1,costhetaB)*legendreP(k2,costhetaC)*legendreP(k3,costhetaA);
        			P =    P + legendreP(k1,costhetaC)*legendreP(k2,costhetaA)*legendreP(k3,costhetaB);
        			P =    P + legendreP(k1,costhetaC)*legendreP(k2,costhetaB)*legendreP(k3,costhetaA);
        			
        			Vexp = Vexp + A[k1][k2][k3]*Math.exp(-alpha[k1][k2][k3]*(RAB+RBC+RAC))*P;
        			
        			//System.out.println(k1 + " " + k2 + " " + k3 + "  " + alpha[k1][k2][k3]);
        		}
        	}
        }
        
        ///////////////////////////////////////////////////////////////////
        //V3Disp
        ///////////////////////////////////////////////////////////////////
        
        double R1=RAB; double R2 = RAC; double R3 = RBC;
        double theta1=thetaA; double theta2=thetaB; double theta3=thetaC;
        
        // W111: all R and theta treated equivalently
        double W111 = 3.0*Math.pow(R1,-3)*Math.pow(R2,-3)*Math.pow(R3,-3);
               W111 = W111*(1.0 + (3.0*Math.cos(theta1)*Math.cos(theta2)*Math.cos(theta3) ));
        double D111 = getD(beta3[1][1][1],R1,3)*getD(beta3[1][1][1],R2,3)*getD(beta3[1][1][1],R3,3);
       
        // W112: R1 and theta3 unique; let each R and theta get to be R1 and theta3
        R1=RAB; R2 = RAC;R3 = RBC;
        theta1=thetaA; theta2=thetaB; theta3=thetaC;
        double W112 = 3.0/16.0*Math.pow(R1,-3)*Math.pow(R2,-4)*Math.pow(R3,-4);
        double term112 = (9.0*Math.cos(theta3)-25.0*Math.cos(3.0*theta3));
               term112 = term112 + (6.0*Math.cos(theta1-theta2)*(3.0+5.0*Math.cos(2.0*theta3)));
               W112 = W112*term112;
        double D112 = getD(beta3[1][1][2],R1,3)*getD(beta3[1][1][2],R2,4)*getD(beta3[1][1][2],R3,4);
        // W112: Let RBC get the chance to be R1 (R2 and R3 treated equivalently)
        R1=RBC; R2 = RAB; R3 = RAC;
   	    theta1=thetaB; theta2=thetaC; theta3=thetaA;
        double W121 = 3.0/16.0*Math.pow(R1,-3)*Math.pow(R2,-4)*Math.pow(R3,-4);
        double term121 = (9.0*Math.cos(theta3)-25.0*Math.cos(3.0*theta3));
               term121 = term121 + (6.0*Math.cos(theta1-theta2)*(3.0+5.0*Math.cos(2.0*theta3)));
               W121 = W121*term121;
        double D121 = getD(beta3[1][1][2],R1,3)*getD(beta3[1][1][2],R2,4)*getD(beta3[1][1][2],R3,4);
        // W112: Let RAC get the chance to be R1 (R2 and R3 treated equivalently)
        R1=RAC; R2 = RAB; R3 = RBC;
        theta1=thetaA; theta2=thetaC; theta3=thetaB;
        double W211 = 3.0/16.0*Math.pow(R1,-3)*Math.pow(R2,-4)*Math.pow(R3,-4);
        double term211 = (9.0*Math.cos(theta3)-25.0*Math.cos(3.0*theta3));
               term211 = term211 + (6.0*Math.cos(theta1-theta2)*(3.0+5.0*Math.cos(2.0*theta3)));
               W211 = W211*term211;
        double D211 = getD(beta3[1][1][2],R1,3)*getD(beta3[1][1][2],R2,4)*getD(beta3[1][1][2],R3,4);
       
        // W122: R2 and theta2 unique; let each R and theta get to be R2 and theta1
        // R2 must be side opposite theta2
        R1=RAC; R2 = RBC; R3 = RAB;
        theta1=thetaC; theta2=thetaA; theta3=thetaB;
        double W122 = 15.0/64.0*Math.pow(R1,-4)*Math.pow(R2,-5)*Math.pow(R3,-4);
        double term122 = 3*(Math.cos(theta2)+5.0*Math.cos(3*theta2));
        	   term122 = term122 +   20.0*Math.cos(theta1-theta3)*(1.0-3.0*Math.cos(2*theta2));
        	   term122 = term122 + 70.0*Math.cos(2.0*(theta1-theta3))*Math.cos(theta2);
        	   W122 = W122*term122;
	    double D122 = getD(beta3[1][2][2],R1,4)*getD(beta3[1][2][2],R2,5)*getD(beta3[1][2][2],R3,4);   
	    // W122: Let RBC get the chance to be R1 (R2 and R3 treated equivalently)
        R1=RAC; R2 = RAB; R3 = RBC;
   	    theta1=thetaA; theta2=thetaC; theta3=thetaB;
   	    double W212 = 15.0/64.0*Math.pow(R1,-4)*Math.pow(R2,-5)*Math.pow(R3,-4);
   	    double term212 = 3*(Math.cos(theta2)+5.0*Math.cos(3*theta2));
     	   term212 = term212 +   20.0*Math.cos(theta1-theta3)*(1.0-3.0*Math.cos(2*theta2));
     	   term212 = term212 + 70.0*Math.cos(2.0*(theta1-theta3))*Math.cos(theta2);
     	   W212 = W212*term212;
	    double D212 = getD(beta3[1][2][2],R1,4)*getD(beta3[1][2][2],R2,5)*getD(beta3[1][2][2],R3,4); 
	    // W112: Let RAC get the chance to be R1 (R2 and R3 treated equivalently)
        R1=RBC; R2 = RAC; R3 = RAB;
        theta1=thetaC; theta2=thetaB; theta3=thetaA;
        double W221 = 15.0/64.0*Math.pow(R1,-4)*Math.pow(R2,-5)*Math.pow(R3,-4);
   	    double term221 = 3*(Math.cos(theta2)+5.0*Math.cos(3*theta2));
     	   term221 = term221 +   20.0*Math.cos(theta1-theta3)*(1.0-3.0*Math.cos(2*theta2));
     	   term221 = term221 + 70.0*Math.cos(2.0*(theta1-theta3))*Math.cos(theta2);
     	   W221 = W221*term221;
	    double D221 = getD(beta3[1][2][2],R1,4)*getD(beta3[1][2][2],R2,5)*getD(beta3[1][2][2],R3,4); 
	    
	    // W222: all R and theta treated equivalently
        double W222 = 15.0/128.0*Math.pow(R1,-5)*Math.pow(R2,-5)*Math.pow(R3,-5);
        double term222 = -27.0 + 220.0*Math.cos(theta1)*Math.cos(theta2)*Math.cos(theta3);
        	   // As in Bukowski and Szalewicz 2001
               //term222 = term222 + 490.0*Math.cos(2.0*theta1)*Math.cos(theta2)*Math.cos(theta3);
               // As in Cencek's et al's code
        	   term222 = term222 +   490.0*Math.cos(2.0*theta1)*Math.cos(2.0*theta2)*Math.cos(2.0*theta3);
               term222 = term222 + 175.0*(Math.cos(2.0*(theta1-theta2))+Math.cos(2.0*(theta2-theta3))+Math.cos(2.0*(theta3-theta1)));
               W222 = W222*term222;
        double D222 = getD(beta3[2][2][2],R1,5)*getD(beta3[2][2][2],R2,5)*getD(beta3[2][2][2],R3,5);
        
        // W113: R1 and theta3 unique; let each R and theta get to be R1 and theta3
        // Permutation 1
        R1=RAB; R2 = RBC; R3 = RAC;
        theta1=thetaB; theta2=thetaA; theta3=thetaC;
        double W113 = 5.0/32.0*Math.pow(R1,-3)*Math.pow(R2,-5)*Math.pow(R3,-5);
        double term113 = 9.0 + 8.0*Math.cos(2.0*theta3) - 49.0*Math.cos(4.0*theta3);
               term113 = term113 + 6.0*Math.cos(theta1-theta2)*(9.0*Math.cos(theta3)+7.0*Math.cos(3.0*theta3));
               W113 = W113*term113;
        double D113 = getD(beta3[1][1][3],R1,3)*getD(beta3[1][1][3],R2,5)*getD(beta3[1][1][3],R3,5);
        // Permutation 2
        R1=RBC; R2 = RAB; R3 = RAC;
   	    theta1=thetaB; theta2=thetaC; theta3=thetaA;
	   	double W131 = 5.0/32.0*Math.pow(R1,-3)*Math.pow(R2,-5)*Math.pow(R3,-5);
	    double term131 = 9.0 + 8.0*Math.cos(2.0*theta3) - 49.0*Math.cos(4.0*theta3);
	           term131 = term131 + 6.0*Math.cos(theta1-theta2)*(9.0*Math.cos(theta3)+7.0*Math.cos(3.0*theta3));
	           W131 = W131*term131;
	    double D131 = getD(beta3[1][1][3],R1,3)*getD(beta3[1][1][3],R2,5)*getD(beta3[1][1][3],R3,5);
	    // Permutation 3
	    R1=RAC; R2 = RAB; R3 = RBC;
	    theta1=thetaA; theta2=thetaC; theta3=thetaB;
	    double W311 = 5.0/32.0*Math.pow(R1,-3)*Math.pow(R2,-5)*Math.pow(R3,-5);
	    double term311 = 9.0 + 8.0*Math.cos(2.0*theta3) - 49.0*Math.cos(4.0*theta3);
	           term311 = term311 + 6.0*Math.cos(theta1-theta2)*(9.0*Math.cos(theta3)+7.0*Math.cos(3.0*theta3));
	           W311 = W311*term311;
	    double D311 = getD(beta3[1][1][3],R1,3)*getD(beta3[1][1][3],R2,5)*getD(beta3[1][1][3],R3,5);
	    
        double V3111 =   D111*W111*Z3[1][1][1];
        double V3112 = (D112*W112 + D121*W121 + D211*W211)*Z3[1][1][2];
        double V3122 = (D122*W122 + D212*W212 + D221*W221)*Z3[1][2][2];
        double V3222 = D222*W222*Z3[2][2][2];
        double V3113 = (D113*W113 + D131*W131 + D311*W311)*Z3[1][1][3];
        double V3disp = V3111 + V3112 + V3122 + V3222 + V3113;
        
        if (verbose) {
	        System.out.println();
	        System.out.println("V3 111 " + V3111*KPerHartree + " K");
	        System.out.println("V3 112 " + V3112*KPerHartree + " K");
	        System.out.println("V3 122 " + V3122*KPerHartree + " K");
	        System.out.println("V3 222 " + V3222*KPerHartree + " K");
	        System.out.println("V3 113 " + V3113*KPerHartree + " K");
        }
		
        
        
        
        ///////////////////////////////////////////////////////////////////
        //V4Disp
        ///////////////////////////////////////////////////////////////////
        
   	    // W1111_220=W660: R3 and theta1 unique; let each R and theta get to be R3 and theta1
   	    // Permutation 1
	   	R1=RAB; R2 = RAC; R3 = RBC;
	    theta1=thetaA; theta2=thetaB; theta3=thetaC;
    	double W660 = 9.0*(1.0+Math.cos(theta1)*Math.cos(theta1))*Math.pow(R1,-6)*Math.pow(R2,-6);
    	double D660 = getD(beta4220[1][1][1][1],R1,6)*getD(beta4220[1][1][1][1],R2,6);
    	// Permutation 2
    	R1=RBC; R2 = RAB; R3 = RAC;
   	    theta1=thetaB; theta2=thetaC; theta3=thetaA;
    	double W066 = 9.0*(1.0+Math.cos(theta1)*Math.cos(theta1))*Math.pow(R1,-6)*Math.pow(R2,-6);
    	double D066 = getD(beta4220[1][1][1][1],R1,6)*getD(beta4220[1][1][1][1],R2,6);
    	// Permutation 3
    	R1=RBC; R2 = RAC; R3 = RAB;
	    theta1=thetaC; theta2=thetaB; theta3=thetaA;
    	double W606 = 9.0*(1.0+Math.cos(theta1)*Math.cos(theta1))*Math.pow(R1,-6)*Math.pow(R2,-6);
        double D606 = getD(beta4220[1][1][1][1],R1,6)*getD(beta4220[1][1][1][1],R2,6);
        
        
        // W1112_211=W734: Nothing treated equivalently
        // W1112_211=W734: Permutation 1
        R1=RAB; R2 = RAC; R3 = RBC;
        theta1=thetaA; theta2=thetaB; theta3=thetaC;
        double W734 = 1.0/32.0*Math.pow(R1,-7)*Math.pow(R2,-4)*Math.pow(R3,-3);
        double term734 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term734 = term734 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta1);
               term734 = term734 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2.0*theta2);
               W734 = W734*term734;  
   	    double D734 = getD(beta4211[1][1][1][2],R1,7)*getD(beta4211[1][1][1][2],R2,4)*getD(beta4211[1][1][1][2],R3,3);
   	    // W1112_211=W734: Permutation 2
        R1=RAB; R2 = RBC; R3 = RAC;
        theta1=thetaB; theta2=thetaA; theta3=thetaC;
        double W743 = 1.0/32.0*Math.pow(R1,-7)*Math.pow(R2,-4)*Math.pow(R3,-3);
        double term743 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term743 = term743 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta1);
               term743 = term743 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2.0*theta2);
               W743 = W743*term743;  
   	    double D743 = getD(beta4211[1][1][1][2],R1,7)*getD(beta4211[1][1][1][2],R2,4)*getD(beta4211[1][1][1][2],R3,3);
   	    // W1112_211=W734: Permutation 3
        R1=RBC; R2 = RAC; R3 = RAB;
   	    theta1=thetaC; theta2=thetaB; theta3=thetaA;
   	    double W374 = 1.0/32.0*Math.pow(R1,-7)*Math.pow(R2,-4)*Math.pow(R3,-3);
        double term374 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term374 = term374 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta1);
               term374 = term374 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2.0*theta2);
               W374 = W374*term374;  
	    double D374 = getD(beta4211[1][1][1][2],R1,7)*getD(beta4211[1][1][1][2],R2,4)*getD(beta4211[1][1][1][2],R3,3);
	    // W1112_211=W734: Permutation 4
        R1=RBC; R2 = RAB; R3 = RAC;
   	    theta1=thetaB; theta2=thetaC; theta3=thetaA;
   	    double W347 = 1.0/32.0*Math.pow(R1,-7)*Math.pow(R2,-4)*Math.pow(R3,-3);
        double term347 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term347 = term347 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta1);
               term347 = term347 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2.0*theta2);
               W347 = W347*term347;  
	    double D347 = getD(beta4211[1][1][1][2],R1,7)*getD(beta4211[1][1][1][2],R2,4)*getD(beta4211[1][1][1][2],R3,3);
	    // W1112_211=W734: Permutation 5
        R1=RAC; R2 = RAB; R3 = RBC;
   	    theta1=thetaA; theta2=thetaC; theta3=thetaB;
   	    double W473 = 1.0/32.0*Math.pow(R1,-7)*Math.pow(R2,-4)*Math.pow(R3,-3);
        double term473 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term473 = term473 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta1);
               term473 = term473 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2.0*theta2);
               W473 = W473*term473;  
	    double D473 = getD(beta4211[1][1][1][2],R1,7)*getD(beta4211[1][1][1][2],R2,4)*getD(beta4211[1][1][1][2],R3,3);
	    // W1112_211=W734: Permutation 6
        R1=RAC; R2 = RBC; R3 = RAB;
   	    theta1=thetaC; theta2=thetaA; theta3=thetaB;
   	    double W437 = 1.0/32.0*Math.pow(R1,-7)*Math.pow(R2,-4)*Math.pow(R3,-3);
        double term437 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term437 = term437 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta1);
               term437 = term437 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2.0*theta2);
               W437 = W437*term437;  
	    double D437 = getD(beta4211[1][1][1][2],R1,7)*getD(beta4211[1][1][1][2],R2,4)*getD(beta4211[1][1][1][2],R3,3);
   	    
   	    // W1121_211=W644, R1 and theta3 are unique
	    // Permutation 1
	    R1=RAB; R2 = RAC; R3 = RBC;
        theta1=thetaA; theta2=thetaB; theta3=thetaC;
   	    double W644 = 1.0/32.0*Math.pow(R1,-6)*Math.pow(R2,-4)*Math.pow(R3,-4);
   	    double term644 = -111.0*Math.cos(theta3) - 750.0*Math.cos(3.0*theta3);
               term644 = term644 + 180.0*Math.cos(theta1+theta2) + 108.0*Math.cos(theta1-theta2);
               term644 = term644 - 90.0*Math.cos(theta3-2.0*theta1) - 90.0*Math.cos(theta3-2.0*theta2);
                  W644 = W644*term644;
	    double D644 = getD(beta4211[1][1][2][1],R1,6)*getD(beta4211[1][1][2][1],R2,4)*getD(beta4211[1][1][2][1],R3,4); 
	    // Permutation 2
	    R1=RAC; R2 = RAB; R3 = RBC;
        theta1=thetaA; theta2=thetaC; theta3=thetaB;
   	    double W464 = 1.0/32.0*Math.pow(R1,-6)*Math.pow(R2,-4)*Math.pow(R3,-4);
   	    double term464 = -111.0*Math.cos(theta3) - 750.0*Math.cos(3.0*theta3);
               term464 = term464 + 180.0*Math.cos(theta1+theta2) + 108.0*Math.cos(theta1-theta2);
               term464 = term464 - 90.0*Math.cos(theta3-2.0*theta1) - 90.0*Math.cos(theta3-2.0*theta2);
                  W464 = W464*term464;
	    double D464 = getD(beta4211[1][1][2][1],R1,6)*getD(beta4211[1][1][2][1],R2,4)*getD(beta4211[1][1][2][1],R3,4); 
	    // Permutation 3
	    R1=RBC; R2 = RAB; R3 = RAC;
        theta1=thetaB; theta2=thetaC; theta3=thetaA;
   	    double W446 = 1.0/32.0*Math.pow(R1,-6)*Math.pow(R2,-4)*Math.pow(R3,-4);
   	    double term446 = -111.0*Math.cos(theta3) - 750.0*Math.cos(3.0*theta3);
               term446 = term446 + 180.0*Math.cos(theta1+theta2) + 108.0*Math.cos(theta1-theta2);
               term446 = term446 - 90.0*Math.cos(theta3-2.0*theta1) - 90.0*Math.cos(theta3-2.0*theta2);
                  W446 = W446*term446;
	    double D446 = getD(beta4211[1][1][2][1],R1,6)*getD(beta4211[1][1][2][1],R2,4)*getD(beta4211[1][1][2][1],R3,4); 
	    
	    // W2111_211=W833, R1 and theta3 are unique
	    // Permutation 1
	    R1=RAB; R2 = RAC; R3 = RBC;
        theta1=thetaA; theta2=thetaB; theta3=thetaC;
   	    double W833 = -9.0/2.0*Math.pow(R1,-8)*Math.pow(R2,-3)*Math.pow(R3,-3);
   	    double term833 = Math.cos(2.0*theta1) + Math.cos(2.0*theta2) + 6.0*Math.cos(2.0*theta3);
                  W833 = W833*term833;
	    double D833 = getD(beta4211[2][1][1][1],R1,8)*getD(beta4211[2][1][1][1],R2,3)*getD(beta4211[2][1][1][1],R3,3);  
	    // Permutation 2
	    R1=RAC; R2 = RAB; R3 = RBC;
        theta1=thetaA; theta2=thetaC; theta3=thetaB;
   	    double W383 = -9.0/2.0*Math.pow(R1,-8)*Math.pow(R2,-3)*Math.pow(R3,-3);
   	    double term383 = Math.cos(2.0*theta1) + Math.cos(2.0*theta2) + 6.0*Math.cos(2.0*theta3);
                  W383 = W383*term383;
	    double D383 = getD(beta4211[2][1][1][1],R1,8)*getD(beta4211[2][1][1][1],R2,3)*getD(beta4211[2][1][1][1],R3,3);
	    // Permutation 3
	    R1=RBC; R2 = RAB; R3 = RAC;
        theta1=thetaB; theta2=thetaC; theta3=thetaA;
   	    double W338 = -9.0/2.0*Math.pow(R1,-8)*Math.pow(R2,-3)*Math.pow(R3,-3);
   	    double term338 = Math.cos(2.0*theta1) + Math.cos(2.0*theta2) + 6.0*Math.cos(2.0*theta3);
                  W338 = W338*term338;
	    double D338 = getD(beta4211[2][1][1][1],R1,8)*getD(beta4211[2][1][1][1],R2,3)*getD(beta4211[2][1][1][1],R3,3);
               
	    // W1211_220=W770, R3 and theta1 are unique
	    // Permutation 1
	    R1=RAB; R2 = RAC; R3 = RBC;
        theta1=thetaA; theta2=thetaB; theta3=thetaC;
	    double W770 = -1.0/64.0*Math.pow(R1,-7)*Math.pow(R2,-7)*(1485.0*Math.cos(theta1)+ 384.0*Math.cos(3.0*theta1));
        double D770 = getD(beta4220[1][2][1][1],R1,7)*getD(beta4220[1][2][1][1],R2,7);
        // Permutation 2
	    R1=RAC; R2 = RBC; R3 = RAB;
        theta1=thetaC; theta2=thetaA; theta3=thetaB;
	    double W707 = -1.0/64.0*Math.pow(R1,-7)*Math.pow(R2,-7)*(1485.0*Math.cos(theta1)+ 384.0*Math.cos(3.0*theta1));
        double D707 = getD(beta4220[1][2][1][1],R1,7)*getD(beta4220[1][2][1][1],R2,7);
        // Permutation 3
	    R1=RBC; R2 = RAB; R3 = RAC;
        theta1=thetaB; theta2=thetaC; theta3=thetaA;
        double W077 = -1.0/64.0*Math.pow(R1,-7)*Math.pow(R2,-7)*(1485.0*Math.cos(theta1)+ 384.0*Math.cos(3.0*theta1));
        double D077 = getD(beta4220[1][2][1][1],R1,7)*getD(beta4220[1][2][1][1],R2,7);
        
       
        // W2111_220=W860: Nothing treated equivalently
        // W2111_220=W860: Permutation 1
        R1=RAB; R2 = RAC; R3 = RBC;
        theta1=thetaA; theta2=thetaB; theta3=thetaC;
        double W860 = 0.25*Math.pow(R1,-8)*Math.pow(R2,-6)*(369.0 + 288.0*Math.cos(theta1)*Math.cos(theta1));
        double D860 = getD(beta4220[2][1][1][1],R1,8)*getD(beta4220[2][1][1][1],R2,6);
        // W2111_220=W860: Permutation 2
        R1=RAB; R2 = RBC; R3 = RAC;
        theta1=thetaB; theta2=thetaA; theta3=thetaC;
        double W806 = 0.25*Math.pow(R1,-8)*Math.pow(R2,-6)*(369.0 + 288.0*Math.cos(theta1)*Math.cos(theta1));
        double D806 = getD(beta4220[2][1][1][1],R1,8)*getD(beta4220[2][1][1][1],R2,6);
        // W2111_220=W860: Permutation 3
        R1=RBC; R2 = RAC; R3 = RAB;
   	    theta1=thetaC; theta2=thetaB; theta3=thetaA;
        double W086 = 0.25*Math.pow(R1,-8)*Math.pow(R2,-6)*(369.0 + 288.0*Math.cos(theta1)*Math.cos(theta1));
        double D086 = getD(beta4220[2][1][1][1],R1,8)*getD(beta4220[2][1][1][1],R2,6);
	    // W2111_220=W860: Permutation 4
        R1=RBC; R2 = RAB; R3 = RAC;
   	    theta1=thetaB; theta2=thetaC; theta3=thetaA;
   	    double W680 = 0.25*Math.pow(R1,-8)*Math.pow(R2,-6)*(369.0 + 288.0*Math.cos(theta1)*Math.cos(theta1));
        double D680 = getD(beta4220[2][1][1][1],R1,8)*getD(beta4220[2][1][1][1],R2,6);
	    // W2111_220=W860: Permutation 5
        R1=RAC; R2 = RAB; R3 = RBC;
   	    theta1=thetaA; theta2=thetaC; theta3=thetaB;
   	    double W608 = 0.25*Math.pow(R1,-8)*Math.pow(R2,-6)*(369.0 + 288.0*Math.cos(theta1)*Math.cos(theta1));
   	    double D608 = getD(beta4220[2][1][1][1],R1,8)*getD(beta4220[2][1][1][1],R2,6);
	    // W2111_220=W860: Permutation 6
        R1=RAC; R2 = RBC; R3 = RAB;
   	    theta1=thetaC; theta2=thetaA; theta3=thetaB;
   	    double W068 = 0.25*Math.pow(R1,-8)*Math.pow(R2,-6)*(369.0 + 288.0*Math.cos(theta1)*Math.cos(theta1));
	    double D068 = getD(beta4220[2][1][1][1],R1,8)*getD(beta4220[2][1][1][1],R2,6);
        
        
        // As in Bukowski and Szalewicz 2001
        /*
        double V4disp =   D1111_211*W1111_211*Z4211[1][1][1][1]; 
        V4disp = V4disp + D1111_220*W1111_220*Z4220[1][1][1][1];

        V4disp = V4disp + D1112_211*W1112_211*Z4211[1][1][1][2];
        
        V4disp = V4disp + D1121_211*W1121_211*Z4211[1][1][2][1];
        
        V4disp = V4disp + D2111_211*W2111_211*Z4211[2][1][1][1];
        V4disp = V4disp + D2111_220*W2111_220*Z4220[2][1][1][1];
        
        V4disp = V4disp + D1211_220*W1211_220*Z4220[1][2][1][1];
        */
        
        //As in Cencek's code
        
        double V4d660 = D660*W660*Z4220[1][1][1][1]; // W1111_220=W660
        V4d660 = V4d660 + D606*W606*Z4220[1][1][1][1];
        V4d660 = V4d660 + D066*W066*Z4220[1][1][1][1];

        double V4d734 = D734*W734*Z4211[1][1][1][2]; // W1112_211=W734
        V4d734 = V4d734 + D743*W743*Z4211[1][1][1][2];
        V4d734 = V4d734 + D347*W347*Z4211[1][1][1][2];
        V4d734 = V4d734 + D374*W374*Z4211[1][1][1][2];
        V4d734 = V4d734 + D473*W473*Z4211[1][1][1][2];
        V4d734 = V4d734 + D437*W437*Z4211[1][1][1][2];
        
        double V4d644 = D644*W644*Z4211[1][1][2][1]; // W1121_211=W644
        V4d644 = V4d644 + D464*W464*Z4211[1][1][2][1];
        V4d644 = V4d644 + D446*W446*Z4211[1][1][2][1];
        
        double V4d833 = D833*W833*Z4211[2][1][1][1]; // W2111_211=W833
        V4d833 = V4d833 + D383*W383*Z4211[2][1][1][1];
        V4d833 = V4d833 + D338*W338*Z4211[2][1][1][1];
        
        double V4d770 = D770*W770*Z4220[1][2][1][1]; // W1211_220=W770
        V4d770 = V4d770 + D707*W707*Z4220[1][2][1][1];
        V4d770 = V4d770 + D077*W077*Z4220[1][2][1][1];
        
        double V4d860 = D860*W860*Z4220[2][1][1][1]; // W2111_220=W860:
        V4d860 = V4d860 + D806*W806*Z4220[2][1][1][1];
        V4d860 = V4d860 + D086*W086*Z4220[2][1][1][1];
        V4d860 = V4d860 + D680*W680*Z4220[2][1][1][1];
        V4d860 = V4d860 + D608*W608*Z4220[2][1][1][1];
        V4d860 = V4d860 + D068*W068*Z4220[2][1][1][1];
        
        double V4disp = V4d660 + V4d734 + V4d644 + V4d833 + V4d770 + V4d860;
        
        if (verbose) {
	        System.out.println();
	        System.out.println("V4d660 " + V4d660*KPerHartree+ " K");
	        System.out.println("V4d734 " + V4d734*KPerHartree+ " K");
	        System.out.println("V4d644 " + V4d644*KPerHartree+ " K");
	        System.out.println("V4d833 " + V4d833*KPerHartree+ " K");
	        System.out.println("V4d770 " + V4d770*KPerHartree+ " K");
	        System.out.println("V4d860 " + V4d860*KPerHartree+ " K");
	        System.out.println();
	        System.out.println("Vexp " + Vexp*KPerHartree + " K");
	        System.out.println("V3disp " + V3disp*KPerHartree+ " K");
	        System.out.println("V4disp " + V4disp*KPerHartree+ " K");
        }
        return (Vexp+V3disp+V4disp)*KPerHartree; //Kelvin
        
    }
    

    
    public double getD(double beta, double RXY, double nXY) {
    	
    	double D = 1.0;	
    	double factorial = 1.0;
    	for (int n=1;n<=nXY;n++) {
    		
    		factorial = factorial*n;
    		
    		D = D + Math.pow(beta*RXY,n)/factorial;
    	}
    	
    	D = 1.0 - (Math.exp(-beta*RXY)*D);

    	return D;
    }
    
    
    
    public double legendreP (int i, double x) {
    	
    	if (i == 0) {
    		return 1.0;
    	} else if (i == 1) {
    		return x;
    	} else if (i == 2) {
    		//1/2(3x^2-1)	
    		return 0.5*(3.0*x*x-1.0);	
    	} else if (i == 3) {
    		//1/2(5x^3-3x)
    		return 0.5*(5.0*x*x*x-3.0*x);
    	} else if (i == 4) {
    		//1/8(35x^4-30x^2+3)	
    		return 1.0/8.0*(35.0*x*x*x*x-30.0*x*x+3.0);	
    	} else if (i == 5) {
    		// 1/8(63x^5-70x^3+15x)
    		return 1.0/8.0*(63.0*x*x*x*x*x - 70.0*x*x*x + 15.0*x);
    	} else {
    		throw new RuntimeException("Cannot do that order of Legendre polynomial");
    	}
    }
    
    public void setAlpha() {
    	alpha[0][0][0]=1.16406382984624;
        alpha[0][0][1]=0.593775469740520;
        alpha[0][0][2]=0.579991620360140;
        alpha[0][0][3]=2.16547015214329;
        alpha[0][0][4]=2.98314825911820;
        alpha[0][1][1]=1.41509085761783;
        alpha[0][1][2]=1.68952244399482;
        alpha[0][1][3]=4.01020983745692;
        alpha[0][1][4]=0.918500573699908;
        alpha[0][2][2]=1.24672879488806;
        alpha[0][2][3]=2.05681119414852;
        alpha[0][2][4]=1.31995140236607;
        alpha[0][3][3]=3.73423077881244;
        alpha[0][3][4]=2.03094608487497;
        alpha[0][4][4]=2.16349588270954;
        alpha[1][1][1]=1.14951445564366;
        alpha[1][1][2]=1.13500643353835;
        alpha[1][1][3]=0.706124946374114;
        alpha[1][1][4]=3.38198341068258;
        alpha[1][2][2]=0.586893129829570;
        alpha[1][2][3]=2.50706679697449;
        alpha[1][2][4]=0.947942951899726;
        alpha[1][3][3]=1.32340552994745;
        alpha[1][3][4]=0.400000000000000;
        alpha[1][4][4]=3.08785463812471;
        alpha[2][2][2]=0.958505100822341;
        alpha[2][2][3]=1.18369957176966;
        alpha[2][2][4]=2.88051393389548;
        alpha[2][3][3]=3.89009793612663;
        alpha[2][3][4]=1.75941163648810;
        alpha[2][4][4]=0.663027178763952;
        alpha[3][3][3]=3.85457885867751;
        alpha[3][3][4]=1.10245630997729;
        alpha[3][4][4]=0.928262289570450;
        alpha[4][4][4]=1.10072817300826;
       
    }
    
    public void setA() {
    	A[0][0][0]=7.33779142427628;
    	A[0][0][1]=0.207562291096732E-02;
    	A[0][0][2]=0.201257721465354E-02;
    	A[0][0][3]=-162280.921808245;
    	A[0][0][4]=140426930.572053;
    	A[0][1][1]=3676.77127295706;
    	A[0][1][2]=40721.4699374272;
    	A[0][1][3]=-145478302663.879;
    	A[0][1][4]=-2.38575533677150;
    	A[0][2][2]=24.9222823036941;
    	A[0][2][3]=-136123.231454042;
    	A[0][2][4]=-2036.85749653572;
    	A[0][3][3]=65196789563.7813;
    	A[0][3][4]=-246369.934348762;
    	A[0][4][4]=337609.745545195;
    	A[1][1][1]=-303.389041542620;
    	A[1][1][2]=419.847355533508;
    	A[1][1][3]=0.428281282658233E-01;
    	A[1][1][4]=-6376463380.36163;
    	A[1][2][2]=0.121898031059658E-01;
    	A[1][2][3]=25447929.5505372;
    	A[1][2][4]=1.81213617023749;
    	A[1][3][3]=-2092.82734569787;
    	A[1][3][4]=0.681066685456400E-05;
    	A[1][4][4]=-146635205.381890;
    	A[2][2][2]=6.60802549511257;
    	A[2][2][3]=-344.693694202977;
    	A[2][2][4]=143064734.239803;
    	A[2][3][3]=-1008432134063.46;
    	A[2][3][4]=-9312.13937942436;
    	A[2][4][4]=-0.386910250300881E-02;
    	A[3][3][3]=711864422051.715;
    	A[3][3][4]=18.0536606525932;
    	A[3][4][4]=0.922564385936391;
    	A[4][4][4]=-4.52285479839575;

    }

    public void setBeta3() {
    	beta3[1][1][1]=0.850816031004730;
    	beta3[1][1][2]=1.03935993289613;
    	beta3[1][1][3]=2.35163790098234;
    	beta3[1][2][2]=20.0000000000000;
    	beta3[2][2][2]=7.74979337816275;
    }
    
    public void setZ3(){
    	Z3[1][1][1]=0.49311;
    	Z3[1][1][2]=0.92372;
    	Z3[1][1][3]=4.1241;
    	Z3[1][2][2]=1.7377;
    	Z3[2][2][2]=3.2839;
    }
    
    
    public void setBeta4211() {	
    	// Error in data file provided by Cencek et al.
/*    	beta4211[1][1][2][2]=2.22023197004267;
    	beta4211[1][2][2][1]=2.33977220590245;
    	beta4211[2][1][1][2]=1.96782469219456;*/
    	
    	beta4211[1][1][1][2]=2.22023197004267;
    	beta4211[1][1][2][1]=2.33977220590245;
    	beta4211[2][1][1][1]=1.96782469219456;
    }
    
    public void setZ4211(){
    	Z4211[1][1][1][2]=-370.838300778413;
    	Z4211[1][1][2][1]=673.766716043939;
    	Z4211[2][1][1][1]=-553.474291722504;
    }
    
	public void setBeta4220(){
		beta4220[1][1][1][1]=1.76277419240966;
		beta4220[1][2][1][1]=2.13546395662687;
		beta4220[2][1][1][1]=0.959706781068175;

	}
	
	public void setZ4220(){
		Z4220[1][1][1][1]=-15.2910806164061;
		Z4220[1][2][1][1]=158.205832955569;
		Z4220[2][1][1][1]=112.479143795999;
	}
    
    
    /**
     * Sets the nominal bond angle (in radians)
     */
    public void setAngle(double newAngle) {
        angle = newAngle;
    }
    
    /**
     * Returns the nominal bond angle (in radians)
     */
    public double getAngle() {
        return angle;
    }
    
    public Dimension getAngleDimension() {
        return Angle.DIMENSION;
    }

    /**
     * Sets the characteristic energy of the potential
     */
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    /**
     * Returns the characteristic energy of the potential
     */
    public double getEpsilon() {
        return epsilon;
    }
    
    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public Vector[] gradient(IAtomList atoms) {
       throw new RuntimeException("Sorry, no gradient available yet");
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    protected final Vector drAB, drAC, drBC;
    protected double angle;
    protected double epsilon;
    private static final long serialVersionUID = 1L;
    protected final Vector[] gradient;
    public static boolean bigAngle;
    protected double[][][] alpha = new double [5][5][5];
    protected double[][][] A = new double [5][5][5];
    protected double[][][] beta3 = new double [3][3][4];
    protected double[][][] Z3 = new double [3][3][4];
    protected double[][][][] beta4220= new double [3][3][3][3];
    protected double[][][][] Z4220= new double [3][3][3][3];
    protected double[][][][] beta4211= new double [3][3][3][3];
    protected double[][][][] Z4211= new double [3][3][3][3];
    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    public boolean verbose = false;
    
    
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        P3CPSNonAdditiveHeOrig potential = new P3CPSNonAdditiveHeOrig(space);
      
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        
        double a; double U; Vector r0; Vector r1; Vector r2;
        
        System.out.println("Test configurations from Table 1 of Cencek et al. (2009)");
        System.out.println();
        System.out.println("Equilateral triangle 1, rij = 4 a0");  
        a = 4.0*AngstromPerBohrRadius;
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = potential.energy(atoms);

        System.out.println("here    : " + U*1000+ " mK");
        System.out.println("paper   : -56277 mK"); 
        System.out.println("he3fci.f: " +(-56.276964668880169*1000)+" mK"); 
        
        System.out.println();
        
        System.out.println("Equilateral triangle 2, rij = 5.6 a0"); 
        a = 5.6*AngstromPerBohrRadius;
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = potential.energy(atoms);

        System.out.println("here    : " + U*1000+ " mK");
        System.out.println("paper   : -88.31 mK"); 
        System.out.println("he3fci.f: " +(-0.88311765396772574E-01*1000)+" mK"); 
        
        System.out.println();
        
        System.out.println("Equilateral triangle 3, rij = 7 a0"); 
        a = 7.0*AngstromPerBohrRadius;
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = potential.energy(atoms);

        System.out.println("here    : " + U*1000+ " mK");
        System.out.println("paper   : 16.06 mK"); 
        System.out.println("he3fci.f: " +(0.16054998594895231E-01*1000)+" mK"); 
        
        System.out.println();
        
        System.out.println("Line 1, r12 = 5.6 a0, r13 = 11.2 a0, r23 = 5.6 a0");
        a = 5.6*AngstromPerBohrRadius;
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{2 * a, 0, 0});
       
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = potential.energy(atoms);

        System.out.println("here    : " + U*1000 + " mK");
        System.out.println("paper   : -18.59 mK"); 
        System.out.println("he3fci.f: " +(-0.18590407572441018E-01*1000)+" mK"); 
        
        System.out.println();
        System.out.println("Additional Tests");
        System.out.println();
        
        System.out.println("r12=3.0a0, r23=5.0a0, r13=4.0a0");
        a = AngstromPerBohrRadius;
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{3 * a, 0, 0});
        r2 = Vector.of(new double[]{0, 4 * a, 0});
             
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = potential.energy(atoms);

        System.out.println("here    : " + U*1000 + " mK");
        System.out.println("he3fci.f: " +(-41.362849038365589*1000)+" mK"); 
        System.out.println();
        
        System.out.println("r12=6.0a0, r23=5.0a0; r13=5.0a0");
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{6 * a, 0, 0});
        r2 = Vector.of(new double[]{3 * a, 4 * a, 0});
          
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = potential.energy(atoms);
        
        System.out.println("here    : " + U*1000 + " mK");
        System.out.println("he3fci.f: " +(-0.34399437417427647*1000)+" mK"); //millikelvin
        
       
    }
    

}


