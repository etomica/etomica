/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMolecular;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Kelvin;

/** 
 * P2 6-point potential for Nitrogen with anisotropic contribution
 * Reference: 1. Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
 *             phases, JCP 112(15) 6745 (2000)
 *            2. Mulder A, Michels P.J. and Schouten J.A. The importance of the anisotropic energy term for
 *            	the structure of the solid phases of nitrogen, JCP 105(8) 3235 (1996)
 * 
 * @author Tai Boon Tan
 */
public class P2NitrogenAnisotropic extends PotentialMolecular {

	public P2NitrogenAnisotropic(Space space, double rC) {
		super(2, space);
		work = space.makeVector();
		shift = space.makeVector();
		com1 = space.makeVector();
		com2 = space.makeVector();
		
		zk = new Vector[2];
		zl = new Vector[2];
		
		for(int i=0; i<2; i++){
			zk[i] = space.makeVector();
			zl[i] = space.makeVector();
		}
		rkl = space.makeVector();
		
		C = new double[5];
		C[0] = Kelvin.UNIT.toSim( 415.73107);  //[K]
		C[1] = Kelvin.UNIT.toSim(-1446.74414); //[KA^-1]
		C[2] = Kelvin.UNIT.toSim( 2480.73711); //[KA^-2]
		C[3] = Kelvin.UNIT.toSim(-2766.5419); //[KA^-3]`
		C[4] = Kelvin.UNIT.toSim( 1574.2809); //[KA^-4]
		
        chargeP1P1 = chargeP1 * chargeP1;
        chargeP1P2 = chargeP1 * chargeP2;
        chargeP2P2 = chargeP2 * chargeP2;
        
        this.rC = rC;
        
	}

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair){
		double sum = 0.0;
		double r2 = 0.0;

		IMolecule nitrogena = pair.getMolecule(0);
		IMolecule nitrogenb = pair.getMolecule(1);
		
		// to compute the midpoint distance between the two
		Vector pos1 = (nitrogena.getChildList().get(1)).getPosition();
		Vector pos2 = (nitrogenb.getChildList().get(1)).getPosition();
		
		com1.E(pos1);
		com2.E(pos2);
		
		Vector diff1 = space.makeVector();
		Vector diff2 = space.makeVector();
		
		diff1.Ev1Mv2(com1, nitrogena.getChildList().get(0).getPosition());
		diff2.Ev1Mv2(com2, nitrogenb.getChildList().get(0).getPosition());
					
		com1.PEa1Tv1(-0.5, diff1); 		
		com2.PEa1Tv1(-0.5, diff2);
		
	    /*
         *  to check for the nearest image
         *  if it is not nearest image, zeroShift will return 0.0
         */
		
		work.Ev1Mv2(com1, com2);
		shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
		shift.PE(work);
	
		final boolean zeroShift = shift.squared() < 0.1; 
		r2 = work.squared();
		
		if (r2 > rC*rC){ 
			return 0.0;
		}
		//if(r2<1.6) return Double.POSITIVE_INFINITY;
		
		/*
		 * for the point/ atomic assignment
		 * refer to SpeciesN2.java class
		 * 
		 */
        Vector Pa1l = nitrogena.getChildList().get(2).getPosition();
        Vector Pa2l = nitrogena.getChildList().get(3).getPosition();
        Vector Pa1r = nitrogena.getChildList().get(4).getPosition();
        Vector Pa2r = nitrogena.getChildList().get(5).getPosition();
        
        Vector Pb1l = nitrogenb.getChildList().get(2).getPosition();
        Vector Pb2l = nitrogenb.getChildList().get(3).getPosition();
        Vector Pb1r = nitrogenb.getChildList().get(4).getPosition();
        Vector Pb2r = nitrogenb.getChildList().get(5).getPosition();
        
        double r2QQ = 0*2.25;
        
        if (zeroShift) {
    		/*
    		 * 'for' loop for 4 pairs van der Waals interaction between the 
    		 * 	non-bonded atoms between the 2 molecules 
    		 */
    		zk[0].Ev1Mv2((nitrogena.getChildList().get(0)).getPosition(),
    				(nitrogena.getChildList().get(1)).getPosition());
        	zk[0].normalize();
        	zk[1].Ea1Tv1(-1.0, zk[0]);
        	
        	zl[0].Ev1Mv2((nitrogenb.getChildList().get(0)).getPosition(),
        			(nitrogenb.getChildList().get(1)).getPosition());
        	zl[0].normalize();
        	zl[1].Ea1Tv1(-1.0, zl[0]);
        	
    		for (int i=0; i<2; i++){
    			Vector dist = (nitrogenb.getChildList().get(i)).getPosition();
    			  			
    			for (int j=0; j<2; j++){
    				
    				rkl.Ev1Mv2(dist, (nitrogena.getChildList().get(j).getPosition()));
    				double distr2 = rkl.squared();
    				rkl.normalize();
    				
    				double r = Math.sqrt(distr2);
    				r -= calcAnisotropicRho(zk[j], zl[i], rkl);
    				    				
    				if(r > R1){            // R > R1
    					sum += URgtR1(r*r);
    				
    				} else if (r < R1 && r >= R0){  // R1 > R >= R0
    					sum += UR1gtRgteqR0(r*r);
    					
    				} else if (r < R0){   	// R < R0
    					sum += URltR0(r*r);
    					
    				}
    			}
    			
    		}
    		        	
        	//Pa1l
            r2 = Pa1l.Mv1Squared(Pb1l);
            sum += chargeP1P1/Math.sqrt(r2);
            r2 = Pa1l.Mv1Squared(Pb2l);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
            r2 = Pa1l.Mv1Squared(Pb2r);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
            r2 = Pa1l.Mv1Squared(Pb1r);
            sum += chargeP1P1/Math.sqrt(r2);
        
            //Pa2l
            r2 = Pa2l.Mv1Squared(Pb1l);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
            r2 = Pa2l.Mv1Squared(Pb2l);
            sum += chargeP2P2/Math.sqrt(r2);
            r2 = Pa2l.Mv1Squared(Pb2r);
            sum += chargeP2P2/Math.sqrt(r2);
            r2 = Pa2l.Mv1Squared(Pb1r);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
       
            //Pa2r
            r2 = Pa2r.Mv1Squared(Pb1l);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
            r2 = Pa2r.Mv1Squared(Pb2l);
            sum += chargeP2P2/Math.sqrt(r2);
            r2 = Pa2r.Mv1Squared(Pb2r);
            sum += chargeP2P2/Math.sqrt(r2);
            r2 = Pa2r.Mv1Squared(Pb1r);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
       
            //Pa1r
            r2 = Pa1r.Mv1Squared(Pb1l);
            sum += chargeP1P1/Math.sqrt(r2);
            r2 = Pa1r.Mv1Squared(Pb2l);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
            r2 = Pa1r.Mv1Squared(Pb2r);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            sum += chargeP1P2/Math.sqrt(r2);
            r2 = Pa1r.Mv1Squared(Pb1r);
            sum += chargeP1P1/Math.sqrt(r2);
                    
        } 
        
        else {
        	
    		/*
    		 * 'for' loop for 4 pairs van der Waals interaction between the 
    		 * 	non-bonded atoms between the 2 molecules 
    		 * 
    		 * 
    		 * 
    		 */
        	
    		zk[0].Ev1Mv2((nitrogena.getChildList().get(0)).getPosition(),
    				(nitrogena.getChildList().get(1)).getPosition());
        	zk[0].normalize();
        	zk[1].Ea1Tv1(-1.0, zk[0]);
        	
        	zl[0].Ev1Mv2((nitrogenb.getChildList().get(0)).getPosition(),
        			(nitrogenb.getChildList().get(1)).getPosition());
        	zl[0].normalize();
        	zl[1].Ea1Tv1(-1.0, zl[0]);
        	
    		for (int i=0; i<2; i++){
    			Vector dist = (nitrogenb.getChildList().get(i)).getPosition();
    			shift.TE(-1.0);
    			shift.PE(dist);
    			
    			for (int j=0; j<2; j++){
    				
    				rkl.Ev1Mv2(shift, (nitrogena.getChildList().get(j).getPosition()));
    				double distr2 = rkl.squared();
    				rkl.normalize();
    				
    				double r = Math.sqrt(distr2);
    				r -= calcAnisotropicRho(zk[j], zl[i], rkl);
    				
    				if(r > R1){            // R > R1
    					sum += URgtR1(r*r);
    				
    				} else if (r < R1 && r >= R0){  // R1 > R >= R0
    					sum += UR1gtRgteqR0(r*r);
    					
    				} else if (r < R0){   	// R < R0
    					sum += URltR0(r*r);
    					
    				}
    			}
    			shift.ME(dist);
    			shift.TE(-1.0);
    			
    		}
    		
        	shift.TE(-1.0);
        	shift.PE(Pb1l);
            r2 = Pa1l.Mv1Squared(shift);
            shift.ME(Pb1l);
            sum += chargeP1P1/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);      
            r2 = Pa2l.Mv1Squared(shift);
            shift.ME(Pb1l);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);
            r2 = Pa2r.Mv1Squared(shift);
            shift.ME(Pb1l);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);
            r2 = Pa1r.Mv1Squared(shift);
            shift.ME(Pb1l);
            sum += chargeP1P1/Math.sqrt(r2);
        	shift.TE(-1.0);
            
            ////////////
        	shift.TE(-1.0);
            shift.PE(Pb2l);
            r2 = Pa1l.Mv1Squared(shift);
            shift.ME(Pb2l);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            r2 = Pa2l.Mv1Squared(shift);
            shift.ME(Pb2l);
            sum += chargeP2P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            r2 = Pa2r.Mv1Squared(shift);
            shift.ME(Pb2l);
            sum += chargeP2P2/Math.sqrt(r2);
        	shift.TE(-1.0);
            

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            r2 = Pa1r.Mv1Squared(shift);
            shift.ME(Pb2l);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);
            
            
            //////////////////////
        	shift.TE(-1.0);
            shift.PE(Pb2r);
            r2 = Pa1l.Mv1Squared(shift);
            shift.ME(Pb2r);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            r2 = Pa2l.Mv1Squared(shift);
            shift.ME(Pb2r);
            sum += chargeP2P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            r2 = Pa2r.Mv1Squared(shift);
            shift.ME(Pb2r);
            sum += chargeP2P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            r2 = Pa1r.Mv1Squared(shift);
            if (r2 < r2QQ) return Double.POSITIVE_INFINITY;
            shift.ME(Pb2r);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);
            
            /////////////
        	shift.TE(-1.0);
            shift.PE(Pb1r);
            r2 = Pa1l.Mv1Squared(shift);
            shift.ME(Pb1r);
            sum += chargeP1P1/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            r2 = Pa2l.Mv1Squared(shift);
            shift.ME(Pb1r);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            r2 = Pa2r.Mv1Squared(shift);
            shift.ME(Pb1r);
            sum += chargeP1P2/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            r2 = Pa1r.Mv1Squared(shift);
            shift.ME(Pb1r);
            sum += chargeP1P1/Math.sqrt(r2);
        	shift.TE(-1.0);
            
        }
      
        return sum;																					        
	}

    private double URgtR1(double r2){
    	return A1*Math.exp(-alpha1*Math.sqrt(r2)) - B1/(r2*r2*r2);
    
    }
    
    private double UR1gtRgteqR0(double r2){
    	double sumU = 0;
    	double r = Math.sqrt(r2);
    	double RdiffpowN = 1; 
    	
    	for (int i=0; i<=4; i++){
    		sumU += C[i]*RdiffpowN;
    		RdiffpowN *= (r-R0);
    		
    	}
    	
    	return sumU - B1/(r2*r2*r2);
    }
    
    private double URltR0(double r2){
    	return A2*Math.exp(-alpha2*Math.sqrt(r2)) - B1/(r2*r2*r2);
    	
    }
    
    private double calcAnisotropicRho(Vector zk, Vector zl, Vector rkl){
    	/*
    	 * rho1 and rho2 are from Refenrence 2.
    	 */
    	double rho1 = 0.065;
    	double rho2 = 0.07;
    	double zkdotrkl = zk.dot(rkl);
    	double zldotrkl = zl.dot(rkl);
    	
    	zk.ME(zl);
    	double term1 = rho1*zk.dot(rkl);
    	double term2 = 0.5*rho2*(3*zkdotrkl*zkdotrkl + 3*zldotrkl*zldotrkl -1);
    	    	
    	return (term1+term2);
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    private static final long serialVersionUID = 1L;
    
    
	protected Boundary boundary;
	protected final double chargeP1 = ConformationNitrogen.Echarge[SpeciesN2.indexP1left];
	protected final double chargeP2 = ConformationNitrogen.Echarge[SpeciesN2.indexP2left];
	protected final double chargeP1P1, chargeP1P2, chargeP2P2;
	
	protected final double A1 = Kelvin.UNIT.toSim(9.261205e7); //[K] unit
	protected final double alpha1 = 4.037; //[A^-1]
	
	protected final double B1 = Kelvin.UNIT.toSim(1.79e5); // [KA^6]
	protected final double A2 = Kelvin.UNIT.toSim(1.47248e7); //[K]
	protected final double alpha2 = 3.48; //[A^-1]
	protected final double R0 = 3.01006875; //[A]
	protected final double R1 = 3.4494569; //[A]
	
	protected double[] C;
	
	protected final Vector work, shift;
	protected final Vector com1, com2;
	protected final Vector[] zk, zl;
	protected final Vector rkl;
	public double rC =1.0;
}
