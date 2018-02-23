/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMolecular;
import etomica.potential.PotentialMolecularSoft;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.Kelvin;

/** 
 * P2 6-point potential for Nitrogen.  
 * Reference paper: Etters RD. et al, Phys. Rev B 33(12) 1986
 *
 * 
 * 
 * @author Tai Boon Tan
 */
public class P2NitrogenAB extends PotentialMolecular implements PotentialMolecularSoft{

	public P2NitrogenAB(Space space, double rC) {
		super(2, space);
		
		gradient = new Vector[2];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();
		
		work = space.makeVector();
		shift = space.makeVector();
		com1 = space.makeVector();
		com2 = space.makeVector();
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
    		
    		for (int i=0; i<2; i++){
    			Vector dist = (nitrogena.getChildList().get(i)).getPosition();
    			
    			for (int j=0; j<2; j++){
    				double distr2 =dist.Mv1Squared((nitrogenb.getChildList().get(j)).getPosition());
    				if(Math.sqrt(distr2) > R1){            // R > R1
    					sum += URgtR1(distr2);
    				
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					sum += UR1gtRgteqR0(distr2);
    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					sum += URltR0(distr2);
    					
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
    		 */
        	
    		for (int i=0; i<2; i++){
    			Vector dist = (nitrogenb.getChildList().get(i)).getPosition();
    			shift.TE(-1.0);
    			shift.PE(dist);
    			
    			for (int j=0; j<2; j++){
    				double distr2 = (nitrogena.getChildList().get(j)).getPosition().Mv1Squared(shift);
    				
    				if(Math.sqrt(distr2) > R1){            // R > R1
    					sum += URgtR1(distr2);
    				
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					sum += UR1gtRgteqR0(distr2);
    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					sum += URltR0(distr2);
    					
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
      
        return sum*0.5;																					        
	}
    
	public double virial(IMoleculeList pair) {
		
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
		boundary.nearestImage(work);

		Vector[] grad = gradient(pair);
		
		return work.dot(grad[0]);
	}


	public Vector[] gradient(IMoleculeList pair) {
		
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
	
		gradient[0].E(0.0);
		gradient[1].E(0.0);
		
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
        
        if (zeroShift) {
    		/*
    		 * 'for' loop for 4 pairs van der Waals interaction between the 
    		 * 	non-bonded atoms between the 2 molecules 
    		 */
    		
    		for (int i=0; i<2; i++){
    			for (int j=0; j<2; j++){
    				work.Ev1Mv2(nitrogena.getChildList().get(i).getPosition(), nitrogenb.getChildList().get(j).getPosition());
    				r2 = work.squared();
    				
    				if(Math.sqrt(r2) > R1){            // R > R1
    					gradient[0].PEa1Tv1(dURgtR1(r2)/r2, work);
    				
    				} else if (Math.sqrt(r2) < R1 && Math.sqrt(r2) >= R0){  // R1 > R >= R0
    					gradient[0].PEa1Tv1(dUR1gtRgteqR0(r2)/r2, work);
    					
    				} else if (Math.sqrt(r2) < R0){   	// R < R0
    					gradient[0].PEa1Tv1(dURltR0(r2)/r2, work);
    					
    				}
    			}
    			
    		}
    		        	
        	//Pa1l
            work.Ev1Mv2(Pa1l, Pb1l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa1l, Pb2l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa1l, Pb2r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa1l, Pb1r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
        
            //Pa2l
            work.Ev1Mv2(Pa2l, Pb1l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa2l, Pb2l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa2l, Pb2r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa2l, Pb1r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
       
            //Pa2r
            work.Ev1Mv2(Pa2r, Pb1l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa2r, Pb2l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa2r, Pb2r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa2r, Pb1r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
       
            //Pa1r
            work.Ev1Mv2(Pa1r, Pb1l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa1r, Pb2l);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa1r, Pb2r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pa1r, Pb1r);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
                    
        } 
        
        else {
        	
    		/*
    		 * 'for' loop for 4 pairs van der Waals interaction between the 
    		 * 	non-bonded atoms between the 2 molecules 
    		 * 
    		 * 
    		 */
        	
    		for (int i=0; i<2; i++){
    			Vector dist = (nitrogenb.getChildList().get(i)).getPosition();
    			shift.TE(-1.0);
    			shift.PE(dist);
    			
    			for (int j=0; j<2; j++){
    				
    				work.Ev1Mv2(nitrogena.getChildList().get(j).getPosition(), shift);
    				r2 = work.squared();
    				
    				if(Math.sqrt(r2) > R1){            // R > R1
    					gradient[0].PEa1Tv1(dURgtR1(r2)/r2, work);
    				
    				} else if (Math.sqrt(r2) < R1 && Math.sqrt(r2) >= R0){  // R1 > R >= R0
    					gradient[0].PEa1Tv1(dUR1gtRgteqR0(r2)/r2, work);
    					
    				} else if (Math.sqrt(r2) < R0){   	// R < R0
    					gradient[0].PEa1Tv1(dURltR0(r2)/r2, work);
    					
    				}
    			}
    			shift.ME(dist);
    			shift.TE(-1.0);
    			
    		}
    		
        	shift.TE(-1.0);
        	shift.PE(Pb1l);
            work.Ev1Mv2(Pa1l, shift);
            r2 = work.squared();
            shift.ME(Pb1l);
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);      
            work.Ev1Mv2(Pa2l, shift);
            r2 = work.squared();
            shift.ME(Pb1l);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);
            work.Ev1Mv2(Pa2r, shift);
            r2 = work.squared();
            shift.ME(Pb1l);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);
            work.Ev1Mv2(Pa1r, shift);
            r2 = work.squared();
            shift.ME(Pb1l);
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
            ////////////
        	shift.TE(-1.0);
            shift.PE(Pb2l);
            work.Ev1Mv2(Pa1l, shift);
            r2 = work.squared();
            shift.ME(Pb2l);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            work.Ev1Mv2(Pa2l, shift);
            r2 = work.squared();
            shift.ME(Pb2l);
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            work.Ev1Mv2(Pa2r, shift);
            r2 = work.squared();
            shift.ME(Pb2l);
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            work.Ev1Mv2(Pa1r, shift);
            r2 = work.squared();
            shift.ME(Pb2l);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
            
            //////////////////////
        	shift.TE(-1.0);
            shift.PE(Pb2r);
            work.Ev1Mv2(Pa1l, shift);
            r2 = work.squared();
            shift.ME(Pb2r);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            work.Ev1Mv2(Pa2l, shift);
            r2 = work.squared();
            shift.ME(Pb2r);
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            work.Ev1Mv2(Pa2r, shift);
            r2 = work.squared();
            shift.ME(Pb2r);
            gradient[0].PEa1Tv1(-chargeP2P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            work.Ev1Mv2(Pa1r, shift);
            r2 = work.squared();
            shift.ME(Pb2r);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
            /////////////
        	shift.TE(-1.0);
            shift.PE(Pb1r);
            work.Ev1Mv2(Pa1l, shift);
            r2 = work.squared();
            shift.ME(Pb1r);
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            work.Ev1Mv2(Pa2l, shift);
            r2 = work.squared();
            shift.ME(Pb1r);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            work.Ev1Mv2(Pa2r, shift);
            r2 = work.squared();
            shift.ME(Pb1r);
            gradient[0].PEa1Tv1(-chargeP1P2/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            work.Ev1Mv2(Pa1r, shift);
            r2 = work.squared();
            shift.ME(Pb1r);
            gradient[0].PEa1Tv1(-chargeP1P1/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
        }
        
        gradient[1].Ea1Tv1(-1, gradient[0]);
      
		return gradient;
	}


	public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
		// TODO Auto-generated method stub
		return null;
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
    
    
    private double dURgtR1(double r2){
    	double r = Math.sqrt(r2);
    	return -(alpha1*r)*A1*Math.exp(-alpha1*r) - (-6)*B1/(r2*r2*r2);
    
    }
    
    private double dUR1gtRgteqR0(double r2){
    	double sumU = 0;
    	double r = Math.sqrt(r2);
    	double RdiffpowN = 1; 
    	
    	for (int i=1; i<=4; i++){
    		sumU += i*C[i]*RdiffpowN;
    		RdiffpowN *= (r-R0);
    		
    	}
    	
    	return r*sumU - (-6)*B1/(r2*r2*r2);
    }
    
    private double dURltR0(double r2){
    	double r = Math.sqrt(r2);
    	return -(alpha2*r)*A2*Math.exp(-alpha2*r) - (-6)*B1/(r2*r2*r2);
    	
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    private static final long serialVersionUID = 1L;
    
    protected Vector[] gradient;
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
	public double rC, r2;

}
