package etomica.potential;

import java.util.Date;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;

/**
 * Simplified version of ab initio, non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009.
 *  
 * @author Kate Shaul
 */
public class P3CPSNonAdditiveHeSimplified extends Potential implements PotentialSoft, IPotentialAtomicMultibody {

    public P3CPSNonAdditiveHeSimplified(ISpace space) {
        super(3, space);
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
        gradient = new IVectorMutable[3];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        gradient[2] = space.makeVector();

    }

    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public void setNullRegionMethod(int nullRegionMethod) {
        this.nullRegionMethod = nullRegionMethod;
    }

    public double energy(IAtomList atomSet) {
        IAtom atomA = atomSet.getAtom(0);
        IAtom atomB = atomSet.getAtom(1);
        IAtom atomC = atomSet.getAtom(2);

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
        
        if (Math.abs(costhetaA) > 1) costhetaA /= Math.abs(costhetaA);
        if (Math.abs(costhetaB) > 1) costhetaB /= Math.abs(costhetaB);
        if (Math.abs(costhetaC) > 1) costhetaC /= Math.abs(costhetaC);

        double Vexp = 0;

        double Rsum = RAB+RBC+RAC;
        
        //for (int k3 = 0; k3<=0; k3++) {
        /*
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
        			//Vexp += P;
        			
        			//System.out.println(k1 + " " + k2 + " " + k3 + "  " + alpha[k1][k2][k3]);
        		}
        	}
        }
        */
        
        double PA = 0.5*(3.0*costhetaA*costhetaA-1.0);	
        double PB = 0.5*(3.0*costhetaB*costhetaB-1.0);	
        double PC = 0.5*(3.0*costhetaC*costhetaC-1.0);	
    			
    	Vexp = 3.735269915914973*Math.exp(-1.465020118414995*Rsum)*6*PA*PB*PC;
        	
        
        
        ///////////////////////////////////////////////////////////////////
        //V3Disp
        ///////////////////////////////////////////////////////////////////

    	double beta = 1.941680323996295;
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
        double V3disp = D_RAB*D_RAC*D_RBC*3.0*(1.0 + (3.0*costhetaA*costhetaB*costhetaC ))/(prod*prod*prod)*0.49311;
        

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

    public IVector[] gradient(IAtomList atoms) {
       throw new RuntimeException("Sorry, no gradient available yet");
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();

        P3CPSNonAdditiveHeSimplified potential = new P3CPSNonAdditiveHeSimplified(space);
      
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        
        double a; double U; IVector r0; IVector r1; IVector r2;
        boolean test = true;
        if (test) {
	        System.out.println("Test configurations from Table 1 of Cencek et al. (2009)");
	        System.out.println();
	        System.out.println("Equilateral triangle 1, rij = 4 a0");  
	        a = BohrRadius.UNIT.toSim(4.0);
	        r0 = space.makeVector(new double[] {0,0,0});
	        r1 = space.makeVector(new double[] {a,0,0});
	        r2 = space.makeVector(new double[] {a/2.0,a/2.0*Math.sqrt(3),0});
        
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
	        r0 = space.makeVector(new double[] {0,0,0});
	        r1 = space.makeVector(new double[] {a,0,0});
	        r2 = space.makeVector(new double[] {a/2.0,a/2.0*Math.sqrt(3),0});
	        
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
	        r0 = space.makeVector(new double[] {0,0,0});
	        r1 = space.makeVector(new double[] {a,0,0});
	        r2 = space.makeVector(new double[] {a/2.0,a/2.0*Math.sqrt(3),0});
	        
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
	        r0 = space.makeVector(new double[] {0,0,0});
	        r1 = space.makeVector(new double[] {a,0,0});
	        r2 = space.makeVector(new double[] {2*a,0,0});
	       
	        
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
	        r0 = space.makeVector(new double[] {0,0,0});
	        r1 = space.makeVector(new double[] {3*a,0,0});
	        r2 = space.makeVector(new double[] {0,4*a,0});
	             
	        atom0.getPosition().E(r0);
	        atom1.getPosition().E(r1);
	        atom2.getPosition().E(r2);

	        U = Kelvin.UNIT.fromSim(potential.energy(atoms));
	
	        System.out.println("simplified    : " + U*1000+ " mK");
	        UCencek = Kelvin.UNIT.fromSim(pCencek.energy(atoms));
	        System.out.println("Cencek  : " + UCencek*1000+ " mK");
	        System.out.println();
	        
	        System.out.println("r12=6.0a0, r23=5.0a0; r13=5.0a0");
	        r0 = space.makeVector(new double[] {0,0,0});
	        r1 = space.makeVector(new double[] {6*a,0,0});
	        r2 = space.makeVector(new double[] {3*a,4*a,0});
	          
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
		            r0 = space.makeVector(new double[] {0,0,0});
		            r1 = space.makeVector(new double[] {a,0,0});
		            r2 = space.makeVector(new double[] {2*a,0,0});
	           
	            
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
        	
        	r0 = space.makeVector(new double[] {0,0,0});
        	atom0.getPosition().E(r0);
        	
        	
        	for (int i=1;i<=5;i++) {

            	double r01 = 1.5*i;
        		r1 = space.makeVector(new double[] {r01,0,0});
        		atom1.getPosition().E(r1);
    		
    	    for (int j=1;j<=5;j++) {
    	    	
    	    	double r02 = 1.5*j;

    	    	for (int k=1;k<=5;k++) {
    		
    	    		 double theta = Math.PI/5*k;
    	    		 double x2 = r02*Math.cos(theta);
    	    		 double y2 = r02*Math.sin(theta);
    	    		 r2 = space.makeVector(new double[] {x2,y2,0});
    	    		 atom2.getPosition().E(r2);
    	    		 
    	    		 double r12 = Math.sqrt((x2-r01)*(x2-r01)+y2*y2);
    	    		 
    	    		 U = Kelvin.UNIT.fromSim(potential.energy(atoms));
    	    		 
    		    	 System.out.println(r01+"  "+r02+"  " +r12+"  "+ U);
    		    	

    	    	}
    				
    	    }
        }
    }
       
    }
    
    protected final IVectorMutable drAB, drAC, drBC;
    protected IBoundary boundary;
    private static final long serialVersionUID = 1L;
    protected final IVectorMutable[] gradient;
    public static boolean bigAngle;
    protected final double[][][] alpha = new double [5][5][5];
    protected final double[][][] A = new double [5][5][5];
    protected final double[][] Rpow = new double[3][9];

    public boolean verbose = false;
    private int nullRegionMethod = 1;
}


