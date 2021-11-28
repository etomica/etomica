/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.BohrRadius;
import etomica.units.Hartree;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Simplified version of ab initio, non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009.
 *  
 * @author Kate Shaul
 */
public class P3CPSNonAdditiveHeSimplified implements Potential3Soft {

    public P3CPSNonAdditiveHeSimplified(Space space) {
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
    }

    public void setNullRegionMethod(int nullRegionMethod) {
        this.nullRegionMethod = nullRegionMethod;
    }
    
    public void setParameters(double temperatureK) {
    	
    	//Energies in Hartrees, and distances in Bohr radii.
    	
    	if (temperatureK < 350){ //use 100 K fit with training set 6
    		A=-371.6663602; 
            alpha=1.35300572; 
            Z=1.00129424; 
            B=-33.79160965; 
            b=1.09675524; 
    	} else { // use 500 K fit with training set 0

    		A=-64.48078148; 
            alpha=1.23015508; 
            Z=1.37752536; 
            B=-9.57192113; 
            b=1.01409157; 
    	}
    
    	
    }
    
    public void setParameters(String file) {
    	 
    	
    	String d = "/usr/users/kate/HeData/potentials/u123NonAddCPS2009Model/5ParameterSimplifications/";
    	//String d = "C:/Users/Kate/Documents/MATLAB/Helium-4/Trimer/";

    	
		int count = 0;
		try{	
			
			FileReader fileReader = new FileReader(d+file);			
			BufferedReader bufReader = new BufferedReader(fileReader);
			String line;
			
			while ((line = bufReader.readLine()) != null) {
				
				params[count] = Double.parseDouble(line); 
				
				count++;
				
			}
			bufReader.close();
		}catch (IOException e){	
			throw new RuntimeException(e);
		}
    
    	fval=params[0];
    	exitflag=params[1];
    	A=params[2]; //Keep as Hartrees
        alpha=params[3]; // Keep as inverse bohr
        Z=params[4]; //Keep as inverse bohr
        B=params[5]; //Keep as Hartrees
        b=params[6]; // Keep as inverse bohr
        
        
        System.out.println("Nonadditive trimer potential: A = "+A+", a = "+alpha+", Z = "+Z+", B = "+B+", b = "+b);
        
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

	public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3) {
		double RAB = Math.sqrt(dr12.squared());
		double RAC = Math.sqrt(dr13.squared());
		double RBC = Math.sqrt(dr23.squared());

		double costhetaA =  dr12.dot(dr13)/(RAB*RAC);
		double costhetaB = -dr12.dot(dr23)/(RAB*RBC);
		double costhetaC =  dr13.dot(dr23)/(RAC*RBC);

		return energy(RAB, RAC, RBC, costhetaA, costhetaB, costhetaC);
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

    public double energy(double[] r2) {
        return u(r2[0], r2[1], r2[2]);
    }
    
    protected double energy(double RAB, double RAC, double RBC, double costhetaA, double costhetaB, double costhetaC) {
        RAB = BohrRadius.UNIT.fromSim(RAB);
        RAC = BohrRadius.UNIT.fromSim(RAC);
        RBC = BohrRadius.UNIT.fromSim(RBC);
        
        if (nullRegionMethod==1) {  //Previously 1 was actually zero...
            if (RAB<3 || RAC<3 || RBC<3) {  
                return 0;
            }
        }
        
        if (nullRegionMethod==2) {
            if (RAB<2.5 || RAC<2.5 || RBC<2.5) {  
                return 0;
            }
        }

        double Vexp = 0;

        double Rsum = RAB+RBC+RAC;
        
    	Vexp = A*Math.exp(-alpha*Rsum)*6.0;
    	Vexp = Vexp + B*Math.exp(-b*Rsum)*6*costhetaA*costhetaB*costhetaC;
        
        ///////////////////////////////////////////////////////////////////
        //V3Disp
        ///////////////////////////////////////////////////////////////////

    	double beta = 0.850816031004730;
    	double D = 1.0;
    	double betaRPowFac = 1.0;
    	for (int n=1;n<=3;n++) {	
    		betaRPowFac *= beta*RAB/n;
    		D = D + betaRPowFac;
    	}
    	double D_RAB = 1.0 - (Math.exp(-beta*RAB)*D);
        	
    	D = 1.0; betaRPowFac = 1.0;
    	for (int n=1;n<=3;n++) {
    		betaRPowFac *= beta*RAC/n;
    		D = D + betaRPowFac;
    	}
        double D_RAC = 1.0 - (Math.exp(-beta*RAC)*D);
        	
    	D = 1.0; betaRPowFac = 1.0;
    	for (int n=1;n<=3;n++) {
    		betaRPowFac *= beta*RBC/n;
    		D = D + betaRPowFac;
    	}
        double D_RBC = 1.0 - (Math.exp(-beta*RBC)*D);
        
        double prod = RAB*RAC*RBC;
        double V3disp = D_RAB*D_RAC*D_RBC*3.0*(1.0 + (3.0*costhetaA*costhetaB*costhetaC ))/(prod*prod*prod)*Z;
    	
    	
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

    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    protected final Vector drAB, drAC, drBC;
    public static boolean bigAngle;

    //Energies in Hartrees, and distances in Bohr radii.
    //Default parameter values fitted for nonadditive classical B3 at 100 K with training set 6 at 100 K
    protected double A=-371.6663602;
    protected double alpha=1.35300572;
    protected double Z=1.00129424; 
    protected double B=-33.79160965;
    protected double b=1.09675524;
    protected double exitflag;
    protected double fval;
    public double[] params = new double[7];//first two are fval and exitflag
    public boolean verbose = false;
    private int nullRegionMethod = 2; // What we have been using so far.
}


