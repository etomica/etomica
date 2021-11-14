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
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

/**
 * Simplified version of ab initio, non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009.
 *  
 * @author Kate Shaul
 */
public class P3CPSNonAdditiveHeSimplified implements PotentialSoft, IPotentialAtomicMultibody, Potential3Soft {

    public P3CPSNonAdditiveHeSimplified(Space space) {
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
        gradient = new Vector[3];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        gradient[2] = space.makeVector();

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

    public Vector[] gradient(IAtomList atoms) {
       throw new RuntimeException("Sorry, no gradient available yet");
    }

    public Vector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        P3CPSNonAdditiveHeSimplified potential = new P3CPSNonAdditiveHeSimplified(space);
        System.out.println(potential.energy(new double[]{4.326850577421106e+00, 4.713606275205238e+03, 4.961477925052809e+03}));
        System.exit(1);
//        potential.setParameters("paramsOriginalSimpler.dat");
        //potential.setParameters("fminconSamplingB3NonAdd/params3Sets_RelErrSqrt_EM18_IG2_1e4Iter_100_K_0.dat");
      
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        
        double a; double U; Vector r0; Vector r1; Vector r2;
        boolean test = true;
        if (test) {
	        
	        System.out.println("Equilateral triangle 1, rij = 4 a0");  
	        a = BohrRadius.UNIT.toSim(4.0);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(new double[]{a, 0, 0});
	        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
        
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        P3CPSNonAdditiveHe pCencek = new P3CPSNonAdditiveHe(space);
	        double UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
	       
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
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
	        
	        System.out.println();
	        
	        System.out.println("Equilateral triangle 2, rij = 10 a0"); 
	        a = BohrRadius.UNIT.toSim(10);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(new double[]{a, 0, 0});
	        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
	        
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        
System.out.println();
	        
	        System.out.println("Equilateral triangle, rij = 2.6 a0"); 
	        a = BohrRadius.UNIT.toSim(2.6);
	        r0 = Vector.of(new double[]{0, 0, 0});
	        r1 = Vector.of(new double[]{a, 0, 0});
	        r2 = Vector.of(new double[]{a / 2.0, a / 2.0 * Math.sqrt(3), 0});
	        
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
		            a = BohrRadius.UNIT.toSim(potential.b);
		            r0 = Vector.of(new double[]{0, 0, 0});
		            r1 = Vector.of(new double[]{a, 0, 0});
		            r2 = Vector.of(new double[]{2 * a, 0, 0});
	           
	            
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
        	P2HePCKLJS p2 = new P2HePCKLJS(space);
        	
        	for (int i=1;i<=5;i++) {

            	double r01 = 1.5*i;
        		r1 = Vector.of(new double[]{r01, 0, 0});
        		atom1.getPosition().E(r1);
    		
    	    for (int j=1;j<=5;j++) {
    	    	
    	    	double r02 = 1.5*j;

    	    	for (int k=1;k<=5;k++) {
    		
    	    		 double theta = Math.PI/5*k;
    	    		 double x2 = r02*Math.cos(theta);
    	    		 double y2 = r02*Math.sin(theta);
    	    		 r2 = Vector.of(new double[]{x2, y2, 0});
    	    		 atom2.getPosition().E(r2);
    	    		 
    	    		 double r12 = Math.sqrt((x2-r01)*(x2-r01)+y2*y2);
    	    		 
    	    		 U = Kelvin.UNIT.fromSim(potential.energy(atoms));
    	    		 
    	    		 double Uadd = p2.energy(atoms);
    	    		 
    		    	 System.out.println(r01+"  "+r02+"  " +r12+"  "+ U+"  "+ Uadd);
    		    	

    	    	}
    				
    	    }
        }
    }
       
    }
    
    protected final Vector drAB, drAC, drBC;
    protected final Vector[] gradient;
    public static boolean bigAngle;
    
    
    
    
    
    //Energies in Hartrees, and distances in Bohr radii.
    //Default parameter values fitted for nonadditive classical B3 at 100 K with training set 6 at 100 K
    protected double A=-371.6663602;
    protected double alpha=1.35300572;
    protected double Z=1.00129424; 
    protected double B=-33.79160965;
    protected double b=1.09675524;
    protected final double[][] Rpow = new double[3][9];
    protected double exitflag;
    protected double fval;
    public double[] params = new double[7];//first two are fval and exitflag
    public boolean verbose = false;
    private int nullRegionMethod = 2; // What we have been using so far.
}


