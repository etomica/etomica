/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

/**
 * Ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009.
 *  
 * @author kate, Andrew Schultz
 */
public class P3CPSNonAdditiveHe implements PotentialSoft, IPotentialAtomicMultibody, Potential3Soft {

    public P3CPSNonAdditiveHe(Space space) {
        this(space, 0);
    }

    public P3CPSNonAdditiveHe(Space space, double sigma) {
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
        gradient = new Vector[3];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        gradient[2] = space.makeVector();
        this.sigma = sigma;
        setA();
        setAlpha();
    }

    public void setNullRegionMethod(int nullRegionMethod) {
        this.nullRegionMethod = nullRegionMethod;
    }

    public double energy(IAtomList atomSet) {
        IAtom atomA = atomSet.get(0);
        IAtom atomB = atomSet.get(1);
        IAtom atomC = atomSet.get(2);

        drAB.Ev1Mv2(atomA.getPosition(),atomB.getPosition());
        drAC.Ev1Mv2(atomA.getPosition(),atomC.getPosition());
        drBC.Ev1Mv2(atomB.getPosition(),atomC.getPosition());

        double RAB = Math.sqrt(drAB.squared());
        double RAC = Math.sqrt(drAC.squared());
        double RBC = Math.sqrt(drBC.squared());

        double costhetaA =  drAB.dot(drAC)/(RAB*RAC);
        double costhetaB = -drAB.dot(drBC)/(RAB*RBC);
        double costhetaC =  drAC.dot(drBC)/(RAC*RBC);

        return energy(RAB, RAC, RBC, costhetaA, costhetaB, costhetaC);
    }

    public double energy(double RAB2, double RAC2, double RBC2) {
        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);
        // this fails for R=0, but we bail in that case anyway (below)
        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);
        return energy(RAB, RAC, RBC, costhetaA, costhetaB, costhetaC);
    }

    public double energy(double[] r2) {
        return energy(r2[0], r2[1], r2[2]);
    }
    
    protected double energy(double RAB, double RAC, double RBC, double costhetaA, double costhetaB, double costhetaC) {
        RAB = BohrRadius.UNIT.fromSim(RAB);
        RAC = BohrRadius.UNIT.fromSim(RAC);
        RBC = BohrRadius.UNIT.fromSim(RBC);
        
        if (nullRegionMethod==1) {
            if (RAB<3 || RAC<3 || RBC<3) {                
                return 0;
            }
        }
        if (nullRegionMethod==2) {
            if (RAB<2.5 || RAC<2.5 || RBC<2.5) {  
                return 0;
            }
        }
        
        if (Math.abs(costhetaA) > 1) costhetaA /= Math.abs(costhetaA);
        if (Math.abs(costhetaB) > 1) costhetaB /= Math.abs(costhetaB);
        if (Math.abs(costhetaC) > 1) costhetaC /= Math.abs(costhetaC);

        double Vexp = 0;

        double Rsum = RAB+RBC+RAC;
        
        for (int k3 = 0; k3<=4; k3++) {
            double P3A = legendreP(k3,costhetaA);
            double P3B = legendreP(k3,costhetaB);
            double P3C = legendreP(k3,costhetaC);
            for (int k2 = 0; k2<=k3; k2++) {
                double P2A = legendreP(k2,costhetaA);
                double P2B = legendreP(k2,costhetaB);
                double P2C = legendreP(k2,costhetaC);
                for (int k1 = 0; k1<=k2; k1++) {
                    double P1A = legendreP(k1,costhetaA);
                    double P1B = legendreP(k1,costhetaB);
                    double P1C = legendreP(k1,costhetaC);
        			
        			double P = P1A*(P2B*P3C + P2C*P3B)
        			         + P1B*(P2A*P3C + P2C*P3A)
        			         + P1C*(P2A*P3B + P2B*P3A);
        			
        			Vexp += A[k1][k2][k3]*Math.exp(-alpha[k1][k2][k3]*Rsum)*P;
        			
        			//System.out.println(k1 + " " + k2 + " " + k3 + "  " + alpha[k1][k2][k3]);
        		}
        	}
        }
        
        ///////////////////////////////////////////////////////////////////
        //V3Disp
        ///////////////////////////////////////////////////////////////////

        // 0=A, 1=B, 2=C
        // to match up distances with angles, take
        // 0=BC, 1=AC, 2=AB
        // so, side i is opposite angle i
        
        
        Rpow[2][0] = RAB;
        Rpow[1][0] = RAC;
        Rpow[0][0] = RBC;
        Rpow[2][1] = 1.0/RAB;
        Rpow[1][1] = 1.0/RAC;
        Rpow[0][1] = 1.0/RBC;
        cosNTheta[0][1] = costhetaA;
        cosNTheta[1][1] = costhetaB;
        cosNTheta[2][1] = costhetaC;
        for (int j=0; j<3; j++) {
            for (int i=2; i<9; i++) {
                Rpow[j][i] = Rpow[j][i-1]* Rpow[j][1];
            }
            
            double ct = cosNTheta[j][1];
            sinTheta[j] = Math.sqrt(1-ct*ct);

            double ct2 = ct*ct;
            // cos(N*thetaj)
            cosNTheta[j][2] = 2*ct2 - 1;
            cosNTheta[j][3] = (4*ct2 - 3)*ct;
            cosNTheta[j][4] = (8*ct2*ct2 - 8*ct2 +1);

            double sintheta2j = 2*cosNTheta[j][1]*sinTheta[j];
            for (int k=0; k<j; k++) {
                double v = ct*cosNTheta[k][1] + sinTheta[j]*sinTheta[k];
                // cos(thetaj-thetak)
                cosThetaDiff[j][k] = v;
                cosThetaDiff[k][j] = v;
                // cos(2*(thetaj-thetak))
                cos2ThetaDiff[j][k] = 2*v*v - 1;
                cos2ThetaDiff[k][j] = 2*v*v - 1;
                // sin(2*thetak)
                double sintheta2k = 2*cosNTheta[k][1]*sinTheta[k];
                // cos(thetaj - 2*thetak)
                cos1m2ThetaDiff[j][k] = cosNTheta[j][1]*cosNTheta[k][2] + sinTheta[j]*sintheta2k;
                // cos(thetak - 2*thetaj)
                cos1m2ThetaDiff[k][j] = cosNTheta[j][2]*cosNTheta[k][1] + sinTheta[k]*sintheta2j;
                // cos(thetaj+thetak)
                cosThetaSum[j][k] = ct*cosNTheta[k][1] - sinTheta[j]*sinTheta[k];
                cosThetaSum[k][j] = cosThetaSum[j][k];
            }
        }

        // W111: all R and theta treated equivalently
        double W111 = 3.0*Rpow[2][3]*Rpow[1][3]*Rpow[0][3];
               W111 = W111*(1.0 + (3.0*costhetaA*costhetaB*costhetaC ));
        double D111 = getD(beta3_111,RAB,3)*getD(beta3_111,RAC,3)*getD(beta3_111,RBC,3);
       
        // W112: 0 is distinguishable, 1 and 2 are interchangeable
        double sum112 = doSum112(Rpow[0], Rpow[2], Rpow[1], cosNTheta[0], cosThetaDiff[1][2])
                      + doSum112(Rpow[1], Rpow[2], Rpow[0], cosNTheta[1], cosThetaDiff[0][2])
                      + doSum112(Rpow[2], Rpow[1], Rpow[0], cosNTheta[2], cosThetaDiff[0][1]);
       
        // W122: 0 is distinguishable, 1 and 2 are interchangeable
        double sum122 = doSum122(Rpow[0], Rpow[1], Rpow[2], cosNTheta[0], cosThetaDiff[1][2], cos2ThetaDiff[1][2])
                      + doSum122(Rpow[1], Rpow[0], Rpow[2], cosNTheta[1], cosThetaDiff[0][2], cos2ThetaDiff[0][2])
                      + doSum122(Rpow[2], Rpow[1], Rpow[0], cosNTheta[2], cosThetaDiff[0][1], cos2ThetaDiff[0][1]);

	    // W222: all R and theta treated equivalently
        double W222 = 15.0/128.0*Rpow[2][5]*Rpow[1][5]*Rpow[0][5];
        double term222 = -27.0 + 220.0*costhetaA*costhetaB*costhetaC;
        	   // As in Bukowski and Szalewicz 2001
               //term222 = term222 + 490.0*Math.cos(2.0*theta1)*Math.cos(theta2)*Math.cos(theta3);
               // As in Cencek's et al's code
        	   term222 = term222 +   490.0*cosNTheta[0][2]*cosNTheta[1][2]*cosNTheta[2][2];
               term222 = term222 + 175.0*(cos2ThetaDiff[0][1] + cos2ThetaDiff[1][2] + cos2ThetaDiff[0][2]);
               W222 = W222*term222;
        double D222 = getD(beta3_222,RAB,5)*getD(beta3_222,RAC,5)*getD(beta3_222,RBC,5);
        
        // W113: 0 is distinguishable, 1 and 2 are interchangeable
        double sum113 = doSum113(Rpow[0], Rpow[1], Rpow[2], cosNTheta[0], cosThetaDiff[1][2])
                      + doSum113(Rpow[1], Rpow[0], Rpow[2], cosNTheta[1], cosThetaDiff[0][2])
                      + doSum113(Rpow[2], Rpow[0], Rpow[1], cosNTheta[2], cosThetaDiff[0][1]);
	    
        double V3111 = D111*W111*Z3_111;
        double V3112 = sum112*Z3_112;
        double V3122 = sum122*Z3_122;
        double V3222 = D222*W222*Z3_222;
        double V3113 = sum113*Z3_113;
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
        
   	    // W1111_220=W660: 0 is distinguishable, 1 and 2 are interchangeable
        double sum660 = doSum660(Rpow[1], Rpow[2], cosNTheta[0][1])
                      + doSum660(Rpow[0], Rpow[2], cosNTheta[1][1])
                      + doSum660(Rpow[0], Rpow[1], cosNTheta[2][1]);
        
        
        // W1112_211=W734: Nothing treated equivalently
        double sum734 = doSum734(Rpow[0], Rpow[1], Rpow[2], cosNTheta[2], cosThetaSum[0][1], cosThetaDiff[0][1], cos1m2ThetaDiff[2][0], cos1m2ThetaDiff[2][1])
                      + doSum734(Rpow[0], Rpow[2], Rpow[1], cosNTheta[1], cosThetaSum[0][2], cosThetaDiff[0][2], cos1m2ThetaDiff[1][0], cos1m2ThetaDiff[1][2])
                      + doSum734(Rpow[1], Rpow[0], Rpow[2], cosNTheta[2], cosThetaSum[0][1], cosThetaDiff[0][1], cos1m2ThetaDiff[2][1], cos1m2ThetaDiff[2][0])
                      + doSum734(Rpow[1], Rpow[2], Rpow[0], cosNTheta[0], cosThetaSum[1][2], cosThetaDiff[1][2], cos1m2ThetaDiff[0][1], cos1m2ThetaDiff[0][2])
                      + doSum734(Rpow[2], Rpow[0], Rpow[1], cosNTheta[1], cosThetaSum[0][2], cosThetaDiff[0][2], cos1m2ThetaDiff[1][2], cos1m2ThetaDiff[1][0])
                      + doSum734(Rpow[2], Rpow[1], Rpow[0], cosNTheta[0], cosThetaSum[1][2], cosThetaDiff[1][2], cos1m2ThetaDiff[0][2], cos1m2ThetaDiff[0][1]);
   	    
        // W1121_211=W644, 1 is distinguishable, 2 and 3 are interchangeable
        double sum644 = doSum644(Rpow[0], Rpow[1], Rpow[2], cosNTheta[0], cosThetaSum[1][2], cosThetaDiff[1][2], cos1m2ThetaDiff[0][1], cos1m2ThetaDiff[0][2])
                      + doSum644(Rpow[1], Rpow[0], Rpow[2], cosNTheta[1], cosThetaSum[0][2], cosThetaDiff[0][2], cos1m2ThetaDiff[1][0], cos1m2ThetaDiff[1][2])
                      + doSum644(Rpow[2], Rpow[0], Rpow[1], cosNTheta[2], cosThetaSum[0][1], cosThetaDiff[0][1], cos1m2ThetaDiff[2][0], cos1m2ThetaDiff[2][1]);	    
	    // W2111_211=W833, 1 is distinguishable, 2 and 3 are interchangeable 
        double sum833 = doSum833(Rpow[0], Rpow[1], Rpow[2], cosNTheta[0][2], cosNTheta[1][2], cosNTheta[2][2])
                      + doSum833(Rpow[1], Rpow[0], Rpow[2], cosNTheta[1][2], cosNTheta[0][2], cosNTheta[2][2])
                      + doSum833(Rpow[2], Rpow[0], Rpow[1], cosNTheta[2][2], cosNTheta[0][2], cosNTheta[1][2]);
	    // W1211_220=W770, 0 is distinguishable, 1 and 2 are interchangeable
        double sum770 = doSum770(Rpow[1], Rpow[2], cosNTheta[0])
                      + doSum770(Rpow[0], Rpow[2], cosNTheta[1])
                      + doSum770(Rpow[0], Rpow[1], cosNTheta[2]);

        // W2111_220=W860: Nothing treated equivalently
        double sum860 = doSum860(Rpow[2], Rpow[1], cosNTheta[0][1])
                      + doSum860(Rpow[1], Rpow[2], cosNTheta[0][1])
                      + doSum860(Rpow[0], Rpow[2], cosNTheta[1][1])
                      + doSum860(Rpow[2], Rpow[0], cosNTheta[1][1])
                      + doSum860(Rpow[0], Rpow[1], cosNTheta[2][1])
                      + doSum860(Rpow[1], Rpow[0], cosNTheta[2][1]);

        
        // As in Bukowski and Szalewicz 2001
        /*
        double V4disp =   D1111_211*W1111_211*Z4211[1][1][1][1]; 
        V4disp = V4disp + D1111_220*W1111_220*Z4220_1111;

        V4disp = V4disp + D1112_211*W1112_211*Z4211_1112;
        
        V4disp = V4disp + D1121_211*W1121_211*Z4211_1121;
        
        V4disp = V4disp + D2111_211*W2111_211*Z4211_2111;
        V4disp = V4disp + D2111_220*W2111_220*Z4220_2111;
        
        V4disp = V4disp + D1211_220*W1211_220*Z4220_1211;
        */
        
        //As in Cencek's code
        
        double V4d660 = sum660*Z4220_1111; // W1111_220=W660

        double V4d734 = sum734*Z4211_1112; // W1112_211=W734
        
        double V4d644 = sum644*Z4211_1121; // W1121_211=W644
        
        double V4d833 = sum833*Z4211_2111; // W2111_211=W833
        
        double V4d770 = sum770*Z4220_1211; // W1211_220=W770
        
        double V4d860 = sum860*Z4220_2111; // W2111_220=W860:
        
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

        double u = Hartree.UNIT.toSim(Vexp+V3disp+V4disp);

        if (nullRegionMethod==0) {
            if (RAB<4 && RAC<4 && RBC<4) {
                if (u < 0) {
                    return 0;
                }
            }
        }
        if (sigma != 0) {
            u *= 1 + 0.01 * sigma * Math.signum(u);
        }
        return u;
    }
    

    protected double doSum112(double[] R1pow, double[] R2pow, double[] R3pow, double[] cosNTheta3, double cosThetaDiff12) {
        double W112 = 3.0/16.0*R1pow[3]*R2pow[4]*R3pow[4];
        double term112 = (9.0*cosNTheta3[1]-25.0*cosNTheta3[3]);
               term112 = term112 + (6.0*cosThetaDiff12*(3.0+5.0*cosNTheta3[2]));
               W112 = W112*term112;
        double D112 = getD(beta3_112,R1pow[0],3)*getD(beta3_112,R2pow[0],4)*getD(beta3_112,R3pow[0],4);
        return W112*D112;
    }
    
    protected double doSum122(double[] R1pow, double[] R2pow, double[] R3pow, double[] cosNTheta2, double cosThetaDiff13, double cos2ThetaDiff13) {
        double W122 = 15.0/64.0*R1pow[5]*R2pow[4]*R3pow[4];
        double term122 = 3*(cosNTheta2[1]+5.0*cosNTheta2[3]);
               term122 = term122 + 20.0*cosThetaDiff13*(1.0-3.0*cosNTheta2[2]);
               term122 = term122 + 70.0*cos2ThetaDiff13*cosNTheta2[1];
               W122 = W122*term122;
        double D122 = getD(beta3_122,R2pow[0],4)*getD(beta3_122,R1pow[0],5)*getD(beta3_122,R3pow[0],4);
        return W122*D122;
    }
    
    protected double doSum113(double[] R1, double[] R2, double[] R3, double[] cosNTheta3, double cosThetaDiff12) {
        double W113 = 5.0/32.0*R1[3]*R2[5]*R3[5];
        double term113 = 9.0 + 8.0*cosNTheta3[2] - 49.0*cosNTheta3[4];
               term113 = term113 + 6.0*cosThetaDiff12*(9.0*cosNTheta3[1]+7.0*cosNTheta3[3]);
               W113 = W113*term113;
        double D113 = getD(beta3_113,R1[0],3)*getD(beta3_113,R2[0],5)*getD(beta3_113,R3[0],5);
        return W113*D113;
    }
    
    protected double doSum660(double[] R1, double[] R2, double cosTheta1) {
        double W660 = 9.0*(1.0+cosTheta1*cosTheta1)*R1[6]*R2[6];
        double D660 = getD(beta4220_1111,R1[0],6)*getD(beta4220_1111,R2[0],6);
        return W660*D660;
    }
        
    protected double doSum734(double[] R1, double[] R2, double[] R3, double[] cosNTheta1, double cosThetaSum23, double cosThetaDiff23, double cos1m2Theta13, double cos1m2Theta12) {
        double W734 = 1.0/32.0*R1[7]*R2[4]*R3[3];
        double term734 = -144.0*cosNTheta1[1] + 36.0*cosThetaSum23;
               term734 = term734 + 216.0*cosThetaDiff23 - 120.0*cosNTheta1[3];
               term734 = term734 - 720.0*cos1m2Theta13 - 72.0*cos1m2Theta12;
               W734 = W734*term734;  
        double D734 = getD(beta4211_1112,R1[0],7)*getD(beta4211_1112,R2[0],4)*getD(beta4211_1112,R3[0],3);
        return W734*D734;
    }

    protected double doSum644(double[] R1, double[] R2, double[] R3, double[] cosNTheta3, double cosThetaSum12, double cosThetaDiff12, double cos1m2ThetaDiff31, double cos1m2ThetaDiff32) {
        double W644 = 1.0/32.0*R1[6]*R2[4]*R3[4];
        double term644 = -111.0*cosNTheta3[1] - 750.0*cosNTheta3[3];
               term644 = term644 + 180.0*cosThetaSum12 + 108.0*cosThetaDiff12;
               term644 = term644 - 90.0*cos1m2ThetaDiff31 - 90.0*cos1m2ThetaDiff32;
                  W644 = W644*term644;
        double D644 = getD(beta4211_1121,R1[0],6)*getD(beta4211_1121,R2[0],4)*getD(beta4211_1121,R3[0],4);
        return W644*D644;
    }

    protected double doSum833(double[] R1, double[] R2, double[] R3, double cos2Theta3, double cos2Theta1, double cos2Theta2) {
        double W833 = -9.0/2.0*R1[8]*R2[3]*R3[3];
        double term833 = cos2Theta1 + cos2Theta2 + 6.0*cos2Theta3;
                  W833 = W833*term833;
        double D833 = getD(beta4211_2111,R1[0],8)*getD(beta4211_2111,R2[0],3)*getD(beta4211_2111,R3[0],3);
        return W833*D833;
    }

    protected double doSum770(double[] R1, double[] R2, double[] cosNTheta1) {
        double W770 = -1.0/64.0*R1[7]*R2[7]*(1485.0*cosNTheta1[1]+ 384.0*cosNTheta1[3]);
        double D770 = getD(beta4220_1211,R1[0],7)*getD(beta4220_1211,R2[0],7);
        return W770*D770;
    }

    protected double doSum860(double[] R1, double[] R2, double cosTheta1) {
        double W860 = 0.25*R1[8]*R2[6]*(369.0 + 288.0*cosTheta1*cosTheta1);
        double D860 = getD(beta4220_2111,R1[0],8)*getD(beta4220_2111,R2[0],6);
        return W860*D860;
    }


    public double getD(double beta, double RXY, int nXY) {
    	double D = 1.0;
    	double betaRPowFac = 1.0;

    	// we spend a lot of time in here.  the return value might be used again
    	// once or twice due to permutations.  The intermediate sum from nXY=j
    	// could be used again for nXY=i, so long as i<j
    	
    	// values could be computed via
    	//             D = D + betaPow[n]/RXY[n];
    	// where betaPow are beta^n/n! and RXY are r^-n.  doing this is unhelpful. 
    	double betaRXY = beta*RXY;
    	for (int n=1;n<=nXY;n++) {
    		
    		betaRPowFac *= betaRXY/n;
    		D = D + betaRPowFac;
    	}
    	
    	D = 1.0 - D*Math.exp(-betaRXY);

    	return D;
    }
    
    public double legendreP (int i, double x) {

        switch (i) {
            case 0:
                return 1.0;
            case 1:
                return x;
            case 2:
                //1/2(3x^2-1)	
                return 0.5*(3.0*x*x-1.0);	
            case 3:
                //1/2(5x^3-3x)
                return 0.5*(5.0*x*x-3.0)*x;
            case 4:
        		//1/8(35x^4-30x^2+3)	
        	    double x2 = x*x;
        	    return 1.0/8.0*(35.0*x2*x2-30.0*x2+3.0);	
            case 5:
                // 1/8(63x^5-70x^3+15x)
                x2 = x*x;
                return 1.0/8.0*(63.0*x2*x2 - 70.0*x2 + 15.0)*x;
            default:
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
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    private Vector[] gradient(IAtomList atoms) {
       throw new RuntimeException("Sorry, no gradient available yet");
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        P3CPSNonAdditiveHe potential = new P3CPSNonAdditiveHe(space);
      
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
        a = BohrRadius.UNIT.toSim(4.0);
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = Kelvin.UNIT.fromSim(potential.energy(atoms));

        System.out.println("here    : " + U*1000+ " mK");
        System.out.println("paper   : -56277 mK"); 
        System.out.println("he3fci.f: " +(-56.276964668880169*1000)+" mK"); 
        
        System.out.println();
        
        System.out.println("Equilateral triangle 2, rij = 5.6 a0"); 
        a = BohrRadius.UNIT.toSim(5.6);
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = Kelvin.UNIT.fromSim(potential.energy(atoms));

        System.out.println("here    : " + U*1000+ " mK");
        System.out.println("paper   : -88.31 mK"); 
        System.out.println("he3fci.f: " +(-0.88311765396772574E-01*1000)+" mK"); 
        
        System.out.println();
        
        System.out.println("Equilateral triangle 3, rij = 7 a0"); 
        a = BohrRadius.UNIT.toSim(7.0);
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = Kelvin.UNIT.fromSim(potential.energy(atoms));

        System.out.println("here    : " + U*1000+ " mK");
        System.out.println("paper   : 16.06 mK"); 
        System.out.println("he3fci.f: " +(0.16054998594895231E-01*1000)+" mK"); 
        
        System.out.println();
        
        System.out.println("Line 1, r12 = 5.6 a0, r13 = 11.2 a0, r23 = 5.6 a0");
        a = BohrRadius.UNIT.toSim(5.6);
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{a, 0, 0});
        r2 = Vector.of(new double[]{2 * a, 0, 0});
       
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = Kelvin.UNIT.fromSim(potential.energy(atoms));

        System.out.println("here    : " + U*1000 + " mK");
        System.out.println("paper   : -18.59 mK"); 
        System.out.println("he3fci.f: " +(-0.18590407572441018E-01*1000)+" mK"); 
        
        System.out.println();
        System.out.println("Additional Tests");
        System.out.println();
        
        System.out.println("r12=3.0a0, r23=5.0a0, r13=4.0a0");
        a = BohrRadius.UNIT.toSim(1.0);
        r0 = Vector.of(new double[]{0, 0, 0});
        r1 = Vector.of(new double[]{3 * a, 0, 0});
        r2 = Vector.of(new double[]{0, 4 * a, 0});
             
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);

        U = Kelvin.UNIT.fromSim(potential.energy(atoms));

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

        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
        
        System.out.println("here    : " + U*1000 + " mK");
        System.out.println("he3fci.f: " +(-0.34399437417427647*1000)+" mK"); //millikelvin
        
        for (int i=0;i<50;i++) {
            double b = 1+i*0.1;
            a = BohrRadius.UNIT.toSim(b);
            r0 = Vector.of(new double[]{0, 0, 0});
            r1 = Vector.of(new double[]{a, 0, 0});
            r2 = Vector.of(new double[]{2 * a, 0, 0});
           
            
            atom0.getPosition().E(r0);
            atom1.getPosition().E(r1);
            atom2.getPosition().E(r2);

            U = potential.energy(atoms);

            System.out.println(b+"  " + U);
        }
       
    }
    
    protected final Vector drAB, drAC, drBC;
    protected final Vector[] gradient;
    public static boolean bigAngle;
    protected final double[][][] alpha = new double [5][5][5];
    protected final double[][][] A = new double [5][5][5];
    protected final double[][][] beta3 = new double [3][3][4];
    protected final double[][][] Z3 = new double [3][3][4];
    protected final double[][][][] beta4220= new double [3][3][3][3];
    protected final double[][][][] Z4220= new double [3][3][3][3];
    protected final double[][][][] beta4211= new double [3][3][3][3];
    protected final double[][][][] Z4211= new double [3][3][3][3];
    protected final double[][] Rpow = new double[3][9];
    protected final double[][] cosNTheta = new double[3][5];
    protected final double[] sinTheta = new double[3];
    protected final double[][] cosThetaDiff = new double[3][3];
    protected final double[][] cos2ThetaDiff = new double[3][3];
    protected final double[][] cos1m2ThetaDiff = new double[3][3];
    protected final double[][] cosThetaSum = new double[3][3];
    protected final static double Z3_111=0.49311, Z3_112=0.92372, Z3_113=4.1241, Z3_122=1.7377, Z3_222=3.2839;
    // Error in data file provided by Cencek et al.
    /*      beta4211_1122=2.22023197004267, beta4211_1221=2.33977220590245, beta4211_1221=1.96782469219456;*/
    protected final static double Z4220_1111=-15.2910806164061, Z4220_1211=158.205832955569, Z4220_2111=112.479143795999;
    protected final static double Z4211_1112=-370.838300778413, Z4211_1121=673.766716043939, Z4211_2111=-553.474291722504;
    protected final static double beta3_111=0.850816031004730, beta3_112=1.03935993289613, beta3_113=2.35163790098234, beta3_122=20.0000000000000, beta3_222=7.74979337816275;
    protected final static double beta4211_1112=2.22023197004267, beta4211_1121=2.33977220590245, beta4211_2111=1.96782469219456;
    protected final static double beta4220_1111=1.76277419240966, beta4220_1211=2.13546395662687, beta4220_2111=0.959706781068175;
    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    public boolean verbose = false;
    private int nullRegionMethod = 2; // What we have been using so far.
    protected final double sigma;
}


