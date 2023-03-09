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

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

/**
 * Ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009.
 *  
 * @author kate, Andrew Schultz
 */
public class P3CPSNonAdditiveHeLessSimplified implements IPotential3 {

    public P3CPSNonAdditiveHeLessSimplified(Space space) {
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
    }

    public void setNullRegionMethod(int nullRegionMethod) {
        this.nullRegionMethod = nullRegionMethod;
    }

    public double u(double RAB2, double RAC2, double RBC2) {
        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);
        // this fails for R=0, but we bail in that case anyway (below)
        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);
        return energy(RAB, RAC, RBC, costhetaA, costhetaB, costhetaC);
    }

    public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3) {
        double RAB = Math.sqrt(dr12.squared());
        double RAC = Math.sqrt(dr13.squared());
        double RBC = Math.sqrt(dr23.squared());

        double costhetaA =  dr12.dot(dr13)/(RAB*RAC);
        double costhetaB = -dr12.dot(dr23)/(RAB*RBC);
        double costhetaC =  dr13.dot(dr23)/(RAC*RBC);

        return energy(RAB, RAC, RBC, costhetaA, costhetaB, costhetaC);
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

    public double energy(double[] r2) {
        double RAB2 = r2[0];
        double RAC2 = r2[1];
        double RBC2 = r2[2];
        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);
        // this fails for R=0, but we bail in that case anyway (below)
        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);
        return energy(RAB, RAC, RBC, costhetaA, costhetaB, costhetaC);
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
        
        for (int k3 = 0; k3<=2; k3++) {
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
       
        double V3disp = D111*W111*Z3_111;
        
        if (verbose) {
	        System.out.println();
	        System.out.println("V3 111 " + V3disp*KPerHartree + " K");

        }
		
        
        
        double u = Hartree.UNIT.toSim(Vexp+V3disp); 

        if (nullRegionMethod==0) {
            if (RAB<4 && RAC<4 && RBC<4) {
                if (u < 0) {
                    return 0;
                }
            }
        }
        return u;
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

    	for (int n=1;n<=nXY;n++) {
    		
    		betaRPowFac *= beta*RXY/n;
    		D = D + betaRPowFac;
    	}
    	
    	D = 1.0 - (Math.exp(-beta*RXY)*D);

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
    
    public void setParameters(String file) {
    	 
    	
    	String d = "/usr/users/kate/HeData/potentials/u123NonAddCPS2009Model/21ParameterSimplifications/";
    	//String d = "C:/Users/Kate/Documents/MATLAB/Helium-4/Trimer/fminconBounded/";

    	
		int count = 0;
		try{	
			
			FileReader fileReader = new FileReader(d+file);			
			BufferedReader bufReader = new BufferedReader(fileReader);
			String line;
			
			while ((line = bufReader.readLine()) != null) {
				
				params[count] = Double.parseDouble(line); 
				
				count++;
				
			}
		}catch (IOException e){	
			throw new RuntimeException(e);
		}
    
    	
    	alpha[0][0][0]=params[0];
        alpha[0][0][1]=params[1];
        alpha[0][0][2]=params[2];
        alpha[0][1][1]=params[3];
        alpha[0][1][2]=params[4];
        alpha[0][2][2]=params[5];
        alpha[1][1][1]=params[6];
        alpha[1][1][2]=params[7];
        alpha[1][2][2]=params[8];
        alpha[2][2][2]=params[9];
    	
    	/* original
    	alpha[0][0][0]=1.16406382984624;
        alpha[0][0][1]=0.593775469740520;
        alpha[0][0][2]=0.579991620360140;
        alpha[0][1][1]=1.41509085761783;
        alpha[0][1][2]=1.68952244399482;
        alpha[0][2][2]=1.24672879488806;
        alpha[1][1][1]=1.14951445564366;
        alpha[1][1][2]=1.13500643353835;
        alpha[1][2][2]=0.586893129829570;
        alpha[2][2][2]=0.958505100822341;
        */
    	
    	/* refit at 100 K using fminsearch
    	alpha[0][0][0]=1.2336;
        alpha[0][0][1]=0.3840;
        alpha[0][0][2]=0.6494;
        alpha[0][1][1]=1.7700;
        alpha[0][1][2]=1.4240;
        alpha[0][2][2]=1.3086;
        alpha[1][1][1]=2.7438;
        alpha[1][1][2]=1.6970;
        alpha[1][2][2]=0.8466;
        alpha[2][2][2]=0.9462;
        */
    	
    	/* refit at 500 K with fminsearch, options=statset('TolFun',1e-18, 'TolX',1e-18,'Display', 'iter');
    	alpha[0][0][0]=1.0972e+003;
        alpha[0][0][1]=560.7603;
        alpha[0][0][2]=0.9428;
        alpha[0][1][1]=2.3080e+003;
        alpha[0][1][2]=1.2879;
        alpha[0][2][2]=3.8465e+003;
        alpha[1][1][1]=3.0013e+003;
        alpha[1][1][2]=1.2956e+003;
        alpha[1][2][2]=1.7443;
        alpha[2][2][2]=2.6320;
        */
        
        /* 100-500K with fmincon
    	alpha[0][0][0]=1.13523320;
    	alpha[0][0][1]=0.59103072;
    	alpha[0][0][2]=0.78114668;
    	alpha[0][1][1]=1.63355680;
    	alpha[0][1][2]=1.59632962;
        
    	alpha[0][2][2]=1.94747983;
    	alpha[1][1][1]=1.09503405;
    	alpha[1][1][2]=1.16908945;
    	alpha[1][2][2]=0.61998581;
    	alpha[2][2][2]=1.09799328;
    	*/
    	
    	/*100K (j=0) with fmincon
    	alpha[0][0][0]=1.18817666;
    	alpha[0][0][1]=0.62256917;
    	alpha[0][0][2]=0.75516464;
    	alpha[0][1][1]=1.57228025;
    	alpha[0][1][2]=1.51827171;
    	
    	alpha[0][2][2]=2.11941566;
    	alpha[1][1][1]=1.08261403;
    	alpha[1][1][2]=1.17057771;
    	alpha[1][2][2]=0.72790815;
    	alpha[2][2][2]=1.36990542;
    	*/
    	
    	/*500K (j=0) with fmincon
    	alpha[0][0][0]=1.12665677;
    	alpha[0][0][1]=0.61119644;
    	alpha[0][0][2]=0.77141147;
    	alpha[0][1][1]=1.45100388;
    	alpha[0][1][2]=1.53756199;
    	
    	alpha[0][2][2]=1.75211873;
    	alpha[1][1][1]=1.08742133;
    	alpha[1][1][2]=1.13371281;
    	alpha[1][2][2]=0.75379057;
    	alpha[2][2][2]=1.12571745;
    	*/
    
    	
    	A[0][0][0]=params[10];
    	A[0][0][1]=params[11];
    	A[0][0][2]=params[12];
    	A[0][1][1]=params[13];
    	A[0][1][2]=params[14];
    	A[0][2][2]=params[15];
    	A[1][1][1]=params[16];
    	A[1][1][2]=params[17];
    	A[1][2][2]=params[18];
    	A[2][2][2]=params[19];
    	
    	/* Not Refit
    	A[0][0][0]=7.33779142427628;
    	A[0][0][1]=0.207562291096732E-02;
    	A[0][0][2]=0.201257721465354E-02;
    	A[0][1][1]=3676.77127295706;
    	A[0][1][2]=40721.4699374272;
    	A[0][2][2]=24.9222823036941;
    	A[1][1][1]=-303.389041542620;
    	A[1][1][2]=419.847355533508;
    	A[1][2][2]=0.121898031059658E-01;
    	A[2][2][2]=6.60802549511257;
    	*/
    	
    	/* Refit at 100K using fminsearch
    	A[0][0][0]=5.5655;
    	A[0][0][1]=3.0844e-005;
    	A[0][0][2]=0.0018;
    	A[0][1][1]=3.2078e+003;
    	A[0][1][2]=2.1799e+004;
    	A[0][2][2]=12.3834;
    	A[1][1][1]=-91.8286;
    	A[1][1][2]=454.5106;
    	A[1][2][2]=0.0083;
    	A[2][2][2]=2.1569;
    	*/
    	
    	/*Refit at 500K using fminsearch
    	A[0][0][0]=-3.6426e+003;
    	A[0][0][1]=5.4611;
    	A[0][0][2]=1.2881;
    	A[0][1][1]=-1.3642e+007;
    	A[0][1][2]=2.5607e+003;
    	A[0][2][2]=-2.0488e+005;
    	A[1][1][1]=5.2797e+005;
    	A[1][1][2]=4.4318e+005;
    	A[1][2][2]=-14.7755;
    	A[2][2][2]=5.7432e+003;
    	*/
    	
    	/* 100-500 K with fmincon
    	A[0][0][0]=4.24874900;
    	A[0][0][1]=0.00193830;
    	A[0][0][2]=0.00141854;
    	A[0][1][1]=3676.26905567;
    	A[0][1][2]=40721.45218181;
       
    	A[0][2][2]=24.62908328;
    	A[1][1][1]=92.97191649;
    	A[1][1][2]=375.62101011;
    	A[1][2][2]=0.00623468;
    	A[2][2][2]=5.89995686;
    	*/
    	
    	/*100K (j=0) with fmincon
    	A[0][0][0]=14.67558285;
    	A[0][0][1]=0.00203153;
    	A[0][0][2]=0.00000000;
    	A[0][1][1]=3676.69010301;
    	A[0][1][2]=40721.53235706;
    	
    	A[0][2][2]=26.37844719;
    	A[1][1][1]=59.29141814;
    	A[1][1][2]=413.06114451;
    	A[1][2][2]=0.01457322;
    	A[2][2][2]=7.24175082;
    	*/
    	
    	/* 500K (j=0) with fmincon
    	A[0][0][0]=9.77427499;
    	A[0][0][1]=0.00299417;
    	A[0][0][2]=0.00159213;
    	A[0][1][1]=3676.05631947;
    	A[0][1][2]=40721.47685457;
    	
    	A[0][2][2]=22.74528542;
    	A[1][1][1]=94.51466843;
    	A[1][1][2]=369.56984667;
    	A[1][2][2]=0.02386503;
    	A[2][2][2]=12.96869781;
    	*/


    
    	beta3_111 = params[20];
    	
    	//Original
    	//beta3_111=0.850816031004730;
    	
        // 100 K using fminsearch
        //protected final static double beta3_111=0.8039; //
        // 500 K using fminsearch
        //protected final static double beta3_111 = 7.3884e+003;
        
        // 100-500 K using fmincon
        //protected final static double beta3_111 = 1.35641661;
        
    	//100K (j=0) with fmincon
        //protected final static double beta3_111 = 1.62580974;
        
        //500K (j=0) with fmincon
        //protected final static double beta3_111 = 1.22544936;
    			
    	
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    private Vector[] gradient(IAtomList atoms) {
       throw new RuntimeException("Sorry, no gradient available yet");
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        P3CPSNonAdditiveHeLessSimplified potential = new P3CPSNonAdditiveHeLessSimplified(space);
      
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        
        double a; double U; Vector r0; Vector r1; Vector r2;
        
        boolean test = false;
        if (test) {
	        System.out.println("Test configurations from Table 1 of Cencek et al. (2009)");
	        System.out.println();
	        System.out.println("Equilateral triangle 1, rij = 4 a0");  
	        a = BohrRadius.UNIT.toSim(4.0);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(a, 0, 0);
	        r2 = Vector.of(a / 2.0, a / 2.0 * Math.sqrt(3), 0);
        
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        P3CPSNonAdditiveHeLessSimplified pCencek = new P3CPSNonAdditiveHeLessSimplified(space);
	        double UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
	       
	        System.out.println();
        
	        System.out.println("Equilateral triangle 2, rij = 5.6 a0"); 
	        a = BohrRadius.UNIT.toSim(5.6);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(a, 0, 0);
	        r2 = Vector.of(a / 2.0, a / 2.0 * Math.sqrt(3), 0);
	        
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
	        
	        System.out.println();
	        
	        System.out.println("Equilateral triangle 3, rij = 7 a0"); 
	        a = BohrRadius.UNIT.toSim(7.0);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(a, 0, 0);
	        r2 = Vector.of(a / 2.0, a / 2.0 * Math.sqrt(3), 0);
	        
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        System.out.println();
	        
	        System.out.println("Line 1, r12 = 5.6 a0, r13 = 11.2 a0, r23 = 5.6 a0");
	        a = BohrRadius.UNIT.toSim(5.6);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(a, 0, 0);
	        r2 = Vector.of(2 * a, 0, 0);
	       
	        
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
	        System.out.println();
	        System.out.println("Additional Tests");
	        System.out.println();
        
	        System.out.println("r12=3.0a0, r23=5.0a0, r13=4.0a0");
	        a = BohrRadius.UNIT.toSim(1.0);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(3 * a, 0, 0);
	        r2 = Vector.of(0, 4 * a, 0);
	             
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        System.out.println();
	        
	        System.out.println("r12=6.0a0, r23=5.0a0; r13=5.0a0");
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(6 * a, 0, 0);
	        r2 = Vector.of(3 * a, 4 * a, 0);
	          
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	        
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
	        Date date = new Date();
		       long t0 = date.getTime();
		       for (int i=0;i<1000000;i++) {
		            double b = 1+i*0.1;
		            a = BohrRadius.UNIT.toSim(b);
		            r0 = Vector.of(new double[]{0, 0, 0});
		            r1 = Vector.of(a, 0, 0);
		            r2 = Vector.of(2 * a, 0, 0);
	           
	            
		            atom0.getPosition().E(r0);
		            atom1.getPosition().E(r1);
		            atom2.getPosition().E(r2);
		
		            U = potential.energy(atoms);
		
		            //System.out.println(b+"  " + U);
		        }
		       Date date2 = new Date();
		       long tf = date2.getTime();
		       long tElapsed = tf-t0;
		       System.out.println(tElapsed+"  milliseconds for 1000000 configurations"+"  "+t0+" "+tf);
        } else {
        	
        	r0 = Vector.of(new double[]{0, 0, 0});
        	atom0.getPosition().E(r0);
        	
        	P2HePCKLJS p2 = new P2HePCKLJS();
        	for (int i=1;i<=5;i++) {

            	double r01 = 1.5*i;
        		r1 = Vector.of(r01, 0, 0);
        		atom1.getPosition().E(r1);
    		
    	    for (int j=1;j<=5;j++) {
    	    	
    	    	double r02 = 1.5*j;

    	    	for (int k=1;k<=5;k++) {
    		
    	    		 double theta = Math.PI/5*k;
    	    		 double x2 = r02*Math.cos(theta);
    	    		 double y2 = r02*Math.sin(theta);
    	    		 r2 = Vector.of(x2, y2, 0);
    	    		 atom2.getPosition().E(r2);
    	    		 
    	    		 double r12 = Math.sqrt((x2-r01)*(x2-r01)+y2*y2);
    	    		 
    	    		 U = Kelvin.UNIT.fromSim(potential.energy(atoms));
    	    		 
    	    		 double Uadd = p2.u(r01*r01)+p2.u(r02*r02)+p2.u(r12*r12);
    	    		 
    		    	 System.out.println(r01+"  "+r02+"  " +r12+"  "+ U+"  "+ Uadd);
    		    	


    	    	}
    				
    	    }
        	}
        }
        	
       
    }
    Integer parameters = 0;
    public  double[] params = new double[21];
    protected final Vector drAB, drAC, drBC;
    public static boolean bigAngle;
    protected final double[][][] alpha = new double [5][5][5];
    protected final double[][][] A = new double [5][5][5];
    protected double beta3_111;
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
    protected final static double Z3_111=0.49311;
    
    //original
    //protected final static double beta3_111=0.850816031004730;
    // 100 K using fminsearch
    //protected final static double beta3_111=0.8039; //
    // 500 K using fminsearch
    //protected final static double beta3_111 = 7.3884e+003;
    
    // 100-500 K using fmincon
    //protected final static double beta3_111 = 1.35641661;
    
	//100K (j=0) with fmincon
    //protected final static double beta3_111 = 1.62580974;
    
    //500K (j=0) with fmincon
    //protected final static double beta3_111 = 1.22544936;
    
    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    public boolean verbose = false;
    private int nullRegionMethod = 2; // What we have been using so far.
}


