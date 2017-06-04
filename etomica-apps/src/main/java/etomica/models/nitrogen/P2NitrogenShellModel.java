/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.nitrogen;

import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.space.Vector;
import etomica.potential.PotentialMolecular;
import etomica.potential.PotentialMolecularSoft;
import etomica.space.Space;
import etomica.space.Tensor;

/** 
 * CAUTION!!!!!! The potential DOES NOT really describe the ALPHA and GAMMA structure at
 * 				 NON-ZERO temperature condition!!
 *       		  (through observation [SimulationGraphic]
 *        and pressure calculation[MeterPressureMolecular])
 * 
 * P2 shell-model potential for Nitrogen.  
 *  Reference: 1. Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
 *              	phases, JCP 112(15) 6745 (2000)
 * 			   2. Jordan P.C. et al, Towards phase transferable potential functions: Methodology
 * 					and application to nitrogen, JCP 103(6) 2272 (1995)
 * @author Tai Boon Tan
 */
public class P2NitrogenShellModel extends PotentialMolecular implements PotentialMolecularSoft{

	public P2NitrogenShellModel(Space space, double rC) {
		super(2, space);
		
		gradient = new Vector[2];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();
				
		work = space.makeVector();
		shift = space.makeVector();
		this.rC = rC;
		
        chargeCC = chargeC * chargeC;
        chargeCN = chargeC * chargeN;
        chargeCP = chargeC * chargeP;
        
        chargeNN = chargeN * chargeN;
        chargeNP = chargeN * chargeP;
        
        chargePP = chargeP * chargeP;        
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
		Vector com1 = (nitrogena.getChildList().getAtom(2)).getPosition();
		Vector com2 = (nitrogenb.getChildList().getAtom(2)).getPosition();
		
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
		 * refer to SpeciesN2ShellModel.java class
		 * 
		 */
        Vector Pan1 = nitrogena.getChildList().getAtom(0).getPosition();
        Vector Pan2 = nitrogena.getChildList().getAtom(1).getPosition();
        Vector Pac  = nitrogena.getChildList().getAtom(2).getPosition();
        Vector Pap1 = nitrogena.getChildList().getAtom(3).getPosition();
        Vector Pap2 = nitrogena.getChildList().getAtom(4).getPosition();

        Vector Pbn1 = nitrogenb.getChildList().getAtom(0).getPosition();
        Vector Pbn2 = nitrogenb.getChildList().getAtom(1).getPosition();
        Vector Pbc  = nitrogenb.getChildList().getAtom(2).getPosition();
        Vector Pbp1 = nitrogenb.getChildList().getAtom(3).getPosition();
        Vector Pbp2 = nitrogenb.getChildList().getAtom(4).getPosition();

        
        double r2QQ = 0*2.25;
        
        if (zeroShift) {
    		/*
    		 * 'for' loop for 3 interaction site interaction between the 
    		 * 	non-bonded atoms between the 2 molecules
    		 *  the interactions are designed such a way that 
    		 *  	a. center-center [2 vs 2]
    		 *  	b. center-(-ve)charge [2 vs 3(or)4]
    		 *  
    		 *   so the interaction are [2,2], [2,3], [2,4], [3,2] and [4,2]
    		 *   there is no interaction between the (+ve)charges 
    		 */
    		
    		for (int i=2; i<5; i++){
    			Vector dist = (nitrogena.getChildList().getAtom(i)).getPosition();
    			
    			for (int j=2; j<5; j++){
    				if (i>2 && j!=2) break;
    				
    				double distr2 = dist.Mv1Squared((nitrogenb.getChildList().getAtom(j)).getPosition());
    				
    				if (i==2 && j==2){
    					sum += calcDisperOverlap(distr2, alpha1, epsilon1, delta1);  						
    					
    				} else {
    					sum += calcDisperOverlap(distr2, alpha2, epsilon2, delta2); 
    					
    				}
    			}
    			
    		}  	
        	//Pac
            r2 = Pac.Mv1Squared(Pbc);
            sum += chargeCC/Math.sqrt(r2);
            r2 = Pac.Mv1Squared(Pbp1);
            sum += chargeCP/Math.sqrt(r2);
            r2 = Pac.Mv1Squared(Pbp2);
            sum += chargeCP/Math.sqrt(r2);
            r2 = Pac.Mv1Squared(Pbn1);
            sum += chargeCN/Math.sqrt(r2);
            r2 = Pac.Mv1Squared(Pbn2);
            sum += chargeCN/Math.sqrt(r2);
        
            //Pap1
            r2 = Pap1.Mv1Squared(Pbc);
            sum += chargeCP/Math.sqrt(r2);
            r2 = Pap1.Mv1Squared(Pbp1);
            sum += chargePP/Math.sqrt(r2);
            r2 = Pap1.Mv1Squared(Pbp2);
            sum += chargePP/Math.sqrt(r2);
            r2 = Pap1.Mv1Squared(Pbn1);
            sum += chargeNP/Math.sqrt(r2);
            r2 = Pap1.Mv1Squared(Pbn2);
            sum += chargeNP/Math.sqrt(r2);
       
            //Pap2
            r2 = Pap2.Mv1Squared(Pbc);
            sum += chargeCP/Math.sqrt(r2);
            r2 = Pap2.Mv1Squared(Pbp1);
            sum += chargePP/Math.sqrt(r2);
            r2 = Pap2.Mv1Squared(Pbp2);
            sum += chargePP/Math.sqrt(r2);
            r2 = Pap2.Mv1Squared(Pbn1);
            sum += chargeNP/Math.sqrt(r2);
            r2 = Pap2.Mv1Squared(Pbn2);
            sum += chargeNP/Math.sqrt(r2);
       
            //Pan1
            r2 = Pan1.Mv1Squared(Pbc);
            sum += chargeCN/Math.sqrt(r2);
            r2 = Pan1.Mv1Squared(Pbp1);
            sum += chargeNP/Math.sqrt(r2);
            r2 = Pan1.Mv1Squared(Pbp2);
            sum += chargeNP/Math.sqrt(r2);
            r2 = Pan1.Mv1Squared(Pbn1);
            sum += chargeNN/Math.sqrt(r2);
            r2 = Pan1.Mv1Squared(Pbn2);
            sum += chargeNN/Math.sqrt(r2);
                    
            //Pan2
            r2 = Pan2.Mv1Squared(Pbc);
            sum += chargeCN/Math.sqrt(r2);
            r2 = Pan2.Mv1Squared(Pbp1);
            sum += chargeNP/Math.sqrt(r2);
            r2 = Pan2.Mv1Squared(Pbp2);
            sum += chargeNP/Math.sqrt(r2);
            r2 = Pan2.Mv1Squared(Pbn1);
            sum += chargeNN/Math.sqrt(r2);
            r2 = Pan2.Mv1Squared(Pbn2);
            sum += chargeNN/Math.sqrt(r2);
            
        } 
        
        else {
        	/*
    		 * 'for' loop for 3 interaction site interaction between the 
    		 * 	non-bonded atoms between the 2 molecules
    		 *  the interactions are designed such a way that 
    		 *  	a. center-center [2 vs 2]
    		 *  	b. center-(-ve)charge [2 vs 3(or)4]
    		 *  
    		 *   so the interaction are [2,2], [2,3], [2,4], [3,2] and [4,2]
    		 *   there is no interaction between the (-ve)charges 
    		 */
        	
    		for (int i=2; i<5; i++){
    			Vector dist = (nitrogenb.getChildList().getAtom(i)).getPosition();
    			shift.TE(-1.0);
    			shift.PE(dist);
    			
    			for (int j=2; j<5; j++){
    				if (i>2 && j!=2) break;
    				
    				double distr2 = (nitrogena.getChildList().getAtom(j)).getPosition().Mv1Squared(shift);
    				
    				if (i==2 && j==2){
    					sum += calcDisperOverlap(distr2, alpha1, epsilon1, delta1);  						
    					
    				} else {
    					sum += calcDisperOverlap(distr2, alpha2, epsilon2, delta2); 
    					
    				}
    			}
    			shift.ME(dist);
    			shift.TE(-1.0);
    			
    		}
    		
    		/////// Pbc
    		
        	shift.TE(-1.0);
        	shift.PE(Pbc);
            r2 = Pac.Mv1Squared(shift);
            shift.ME(Pbc);
            sum += chargeCC/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbc);      
            r2 = Pap1.Mv1Squared(shift);
            shift.ME(Pbc);
            sum += chargeCP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbc);
            r2 = Pap2.Mv1Squared(shift);
            shift.ME(Pbc);
            sum += chargeCP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbc);
            r2 = Pan1.Mv1Squared(shift);
            shift.ME(Pbc);
            sum += chargeCN/Math.sqrt(r2);
        	shift.TE(-1.0);
        	
        	shift.TE(-1.0);
            shift.PE(Pbc);
            r2 = Pan2.Mv1Squared(shift);
            shift.ME(Pbc);
            sum += chargeCN/Math.sqrt(r2);
        	shift.TE(-1.0);
            
            ////////// Pbp1
        	shift.TE(-1.0);
            shift.PE(Pbp1);
            r2 = Pac.Mv1Squared(shift);
            shift.ME(Pbp1);
            sum += chargeCP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp1);
            r2 = Pap1.Mv1Squared(shift);
            shift.ME(Pbp1);
            sum += chargePP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp1);
            r2 = Pap2.Mv1Squared(shift);
            shift.ME(Pbp1);
            sum += chargePP/Math.sqrt(r2);
        	shift.TE(-1.0);
            
        	shift.TE(-1.0);
            shift.PE(Pbp1);
            r2 = Pan1.Mv1Squared(shift);
            shift.ME(Pbp1);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);
            
        	shift.TE(-1.0);
            shift.PE(Pbp1);
            r2 = Pan2.Mv1Squared(shift);
            shift.ME(Pbp1);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);
        	
            /////////// Pbp2
        	shift.TE(-1.0);
            shift.PE(Pbp2);
            r2 = Pac.Mv1Squared(shift);
            shift.ME(Pbp2);
            sum += chargeCP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp2);
            r2 = Pap1.Mv1Squared(shift);
            shift.ME(Pbp2);
            sum += chargePP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp2);
            r2 = Pap2.Mv1Squared(shift);
            shift.ME(Pbp2);
            sum += chargePP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp2);
            r2 = Pan1.Mv1Squared(shift);
            shift.ME(Pbp2);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);
            
        	shift.TE(-1.0);
            shift.PE(Pbp2);
            r2 = Pan2.Mv1Squared(shift);
            shift.ME(Pbp2);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);
        	
            /////////// Pbn1
        	shift.TE(-1.0);
            shift.PE(Pbn1);
            r2 = Pac.Mv1Squared(shift);
            shift.ME(Pbn1);
            sum += chargeCN/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn1);
            r2 = Pap1.Mv1Squared(shift);
            shift.ME(Pbn1);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn1);
            r2 = Pap2.Mv1Squared(shift);
            shift.ME(Pbn1);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn1);
            r2 = Pan1.Mv1Squared(shift);
            shift.ME(Pbn1);
            sum += chargeNN/Math.sqrt(r2);
        	shift.TE(-1.0);
        	
        	shift.TE(-1.0);
            shift.PE(Pbn1);
            r2 = Pan2.Mv1Squared(shift);
            shift.ME(Pbn1);
            sum += chargeNN/Math.sqrt(r2);
        	shift.TE(-1.0);
        	
            /////////// Pbn2
        	shift.TE(-1.0);
            shift.PE(Pbn2);
            r2 = Pac.Mv1Squared(shift);
            shift.ME(Pbn2);
            sum += chargeCN/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn2);
            r2 = Pap1.Mv1Squared(shift);
            shift.ME(Pbn2);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn2);
            r2 = Pap2.Mv1Squared(shift);
            shift.ME(Pbn2);
            sum += chargeNP/Math.sqrt(r2);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn2);
            r2 = Pan1.Mv1Squared(shift);
            shift.ME(Pbn2);
            sum += chargeNN/Math.sqrt(r2);
        	shift.TE(-1.0);
        	
        	shift.TE(-1.0);
            shift.PE(Pbn2);
            r2 = Pan2.Mv1Squared(shift);
            shift.ME(Pbn2);
            sum += chargeNN/Math.sqrt(r2);
        	shift.TE(-1.0);
            
        }
      
        return sum;																					        
	}
   
    /*
     * 
     */
	public double virial(IMoleculeList pair) {
		
		IMolecule nitrogena = pair.getMolecule(0);
		IMolecule nitrogenb = pair.getMolecule(1);
		
		// to compute the midpoint distance between the two
		Vector com1 = (nitrogena.getChildList().getAtom(2)).getPosition();
		Vector com2 = (nitrogenb.getChildList().getAtom(2)).getPosition();
		
		work.Ev1Mv2(com1, com2);
		boundary.nearestImage(work);
		
		Vector[] grad = gradient(pair);
		//System.out.println("work: " + work.toString());
		//System.out.println("grad[0]: " + grad[0].toString());
		//System.out.println("grad[1]: " + grad[1].toString());
		//System.exit(1);
		return work.dot(grad[0]);
	}

	/*
	 * Determine the gradient of the molecule
	 * 
	 * It is the energy change with respect to the change in molecular position
	 *  - the change in energy is the summation of the derivative of energy to 
	 *    the atomic position 
	 *  
	 */
	public Vector[] gradient(IMoleculeList pair) {
		
		IMolecule nitrogena = pair.getMolecule(0);
		IMolecule nitrogenb = pair.getMolecule(1);
		
		// to compute the midpoint distance between the two
		Vector com1 = (nitrogena.getChildList().getAtom(2)).getPosition();
		Vector com2 = (nitrogenb.getChildList().getAtom(2)).getPosition();
		
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
		
		//Initial Gradient = 0.0
		gradient[0].E(0.0);
		gradient[1].E(0.0);
		
		//if(r2<1.6) return Double.POSITIVE_INFINITY;
		
		/*
		 * for the point/ atomic assignment
		 * refer to SpeciesN2ShellModel.java class
		 * 
		 */
        Vector Pan1 = nitrogena.getChildList().getAtom(0).getPosition();
        Vector Pan2 = nitrogena.getChildList().getAtom(1).getPosition();
        Vector Pac  = nitrogena.getChildList().getAtom(2).getPosition();
        Vector Pap1 = nitrogena.getChildList().getAtom(3).getPosition();
        Vector Pap2 = nitrogena.getChildList().getAtom(4).getPosition();

        Vector Pbn1 = nitrogenb.getChildList().getAtom(0).getPosition();
        Vector Pbn2 = nitrogenb.getChildList().getAtom(1).getPosition();
        Vector Pbc  = nitrogenb.getChildList().getAtom(2).getPosition();
        Vector Pbp1 = nitrogenb.getChildList().getAtom(3).getPosition();
        Vector Pbp2 = nitrogenb.getChildList().getAtom(4).getPosition();

        if (zeroShift) {
    		/*
    		 * 'for' loop for 3 interaction site interaction between the 
    		 * 	non-bonded atoms between the 2 molecules
    		 *  the interactions are designed such a way that 
    		 *  	a. center-center [2 vs 2]
    		 *  	b. center-(-ve)charge [2 vs 3(or)4]
    		 *  
    		 *   so the interaction are [2,2], [2,3], [2,4], [3,2] and [4,2]
    		 *   there is no interaction between the (+ve)charges 
    		 */
    		
    		for (int i=2; i<5; i++){
    		
    			for (int j=2; j<5; j++){
    				if (i>2 && j!=2) break;
    				
    				work.Ev1Mv2(nitrogena.getChildList().getAtom(i).getPosition(), nitrogenb.getChildList().getAtom(j).getPosition());
    				r2 = work.squared();

    				if (i==2 && j==2){
    					gradient[0].PEa1Tv1(du(r2, alpha1, epsilon1, delta1)/r2, work);   						
    					
    				} else {
    					gradient[0].PEa1Tv1(du(r2, alpha2, epsilon2, delta2)/r2, work);   						
    					
    				}
    			}
    			
    		}  	
        	//Pac
    		work.Ev1Mv2(Pac, Pbc);
    		r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCC/(r2*Math.sqrt(r2)), work);

            work.Ev1Mv2(Pac,Pbp1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
          
            work.Ev1Mv2(Pac, Pbp2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pac,Pbn1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pac,Pbn2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
        
            //Pap1
            work.Ev1Mv2(Pap1, Pbc);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pap1, Pbp1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pap1, Pbp2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pap1, Pbn1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pap1, Pbn2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
       
            //Pap2
            work.Ev1Mv2(Pap2, Pbc);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pap2, Pbp1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
            
            
            work.Ev1Mv2(Pap2, Pbp2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
            
            
            work.Ev1Mv2(Pap2, Pbn1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
            
            
            work.Ev1Mv2(Pap2, Pbn2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
       
            //Pan1
            work.Ev1Mv2(Pan1, Pbc);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan1, Pbp1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan1, Pbp2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan1, Pbn1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan1, Pbn2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
                    
            //Pan2
            work.Ev1Mv2(Pan2, Pbc);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan2, Pbp1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan2, Pbp2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan2, Pbn1);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
            
            work.Ev1Mv2(Pan2, Pbn2);
            r2 = work.squared();
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
                        
        } 
        
        else {
        	/*
    		 * 'for' loop for 3 interaction site interaction between the 
    		 * 	non-bonded atoms between the 2 molecules
    		 *  the interactions are designed such a way that 
    		 *  	a. center-center [2 vs 2]
    		 *  	b. center-(-ve)charge [2 vs 3(or)4]
    		 *  
    		 *   so the interaction are [2,2], [2,3], [2,4], [3,2] and [4,2]
    		 *   there is no interaction between the (-ve)charges 
    		 */
        	
    		for (int i=2; i<5; i++){
    			Vector dist = (nitrogenb.getChildList().getAtom(i)).getPosition();
    			shift.TE(-1.0);
    			shift.PE(dist);
    			
    			for (int j=2; j<5; j++){
    				if (i>2 && j!=2) break;
    				
    				work.Ev1Mv2((nitrogena.getChildList().getAtom(j)).getPosition(), shift);
    				r2 = work.squared();
    				
    				if (i==2 && j==2){
    					gradient[0].PEa1Tv1(du(r2, alpha1, epsilon1, delta1)/r2, work);   						
    					
    				} else {
    					gradient[0].PEa1Tv1(du(r2, alpha2, epsilon2, delta2)/r2, work);   						
    					
    				}
    			}
    			shift.ME(dist);
    			shift.TE(-1.0);
    			
    		}
    		
    		/////// Pbc
    		
        	shift.TE(-1.0);
        	shift.PE(Pbc);
            work.Ev1Mv2(Pac, shift);
            r2 = work.squared();
            shift.ME(Pbc);
            gradient[0].PEa1Tv1(-chargeCC/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbc);      
            work.Ev1Mv2(Pap1, shift);
            r2 = work.squared();
            shift.ME(Pbc);
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbc);
            work.Ev1Mv2(Pap2, shift);
            r2 = work.squared();
            shift.ME(Pbc);
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbc);
            work.Ev1Mv2(Pan1, shift);
            r2 = work.squared();
            shift.ME(Pbc);
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
        	
        	shift.TE(-1.0);
            shift.PE(Pbc);
            work.Ev1Mv2(Pan2, shift);
            r2 = work.squared();
            shift.ME(Pbc);
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
            ////////// Pbp1
        	shift.TE(-1.0);
            shift.PE(Pbp1);
            work.Ev1Mv2(Pac, shift);
            r2 = work.squared();
            shift.ME(Pbp1);
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp1);
            work.Ev1Mv2(Pap1, shift);
            r2 = work.squared();
            shift.ME(Pbp1);
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp1);
            work.Ev1Mv2(Pap2, shift);
            r2 = work.squared();
            shift.ME(Pbp1);
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
        	shift.TE(-1.0);
            shift.PE(Pbp1);
            work.Ev1Mv2(Pan1, shift);
            r2 = work.squared();
            shift.ME(Pbp1);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
        	shift.TE(-1.0);
            shift.PE(Pbp1);
            work.Ev1Mv2(Pan2, shift);
            r2 = work.squared();
            shift.ME(Pbp1);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
        	
            /////////// Pbp2
        	shift.TE(-1.0);
            shift.PE(Pbp2);
            work.Ev1Mv2(Pac, shift);
            r2 = work.squared();
            shift.ME(Pbp2);
            gradient[0].PEa1Tv1(-chargeCP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp2);
            work.Ev1Mv2(Pap1, shift);            
            r2 = work.squared();
            shift.ME(Pbp2);
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp2);
            work.Ev1Mv2(Pap2, shift);            
            r2 = work.squared();
            shift.ME(Pbp2);
            gradient[0].PEa1Tv1(-chargePP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbp2);
            work.Ev1Mv2(Pan1, shift);
            r2 = work.squared();
            shift.ME(Pbp2);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
        	shift.TE(-1.0);
            shift.PE(Pbp2);
            work.Ev1Mv2(Pan2, shift);
            r2 = work.squared();
            shift.ME(Pbp2);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
        	
            /////////// Pbn1
        	shift.TE(-1.0);
            shift.PE(Pbn1);
            work.Ev1Mv2(Pac, shift);
            r2 = work.squared();
            shift.ME(Pbn1);
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn1);
            work.Ev1Mv2(Pap1, shift);
            r2 = work.squared();
            shift.ME(Pbn1);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn1);
            work.Ev1Mv2(Pap2, shift);
            r2 = work.squared();
            shift.ME(Pbn1);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn1);
            work.Ev1Mv2(Pan1, shift);
            r2 = work.squared();
            shift.ME(Pbn1);
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
        	
        	shift.TE(-1.0);
            shift.PE(Pbn1);
            work.Ev1Mv2(Pan2, shift);
            r2 = work.squared();
            shift.ME(Pbn1);
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
        	
            /////////// Pbn2
        	shift.TE(-1.0);
            shift.PE(Pbn2);
            work.Ev1Mv2(Pac, shift);
            r2 = work.squared();
            shift.ME(Pbn2);
            gradient[0].PEa1Tv1(-chargeCN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn2);
            work.Ev1Mv2(Pap1, shift);
            r2 = work.squared();
            shift.ME(Pbn2);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn2);
            work.Ev1Mv2(Pap2, shift);
            r2 = work.squared();
            shift.ME(Pbn2);
            gradient[0].PEa1Tv1(-chargeNP/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pbn2);
            work.Ev1Mv2(Pan1, shift);
            r2 = work.squared();
            shift.ME(Pbn2);
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
        	
        	shift.TE(-1.0);
            shift.PE(Pbn2);
            work.Ev1Mv2(Pan2, shift);
            r2 = work.squared();
            shift.ME(Pbn2);
            gradient[0].PEa1Tv1(-chargeNN/(r2*Math.sqrt(r2)), work);
        	shift.TE(-1.0);
            
        }
      
        gradient[1].Ea1Tv1(-1, gradient[0]);
        
		return gradient;
	}

	public Vector[] gradient(IMoleculeList pair, Tensor pressureTensor) {
		return null;
	}
    
    private double calcDisperOverlap(double r2, double alpha, double epsilon, double delta){
    	double r = Math.sqrt(r2);
    	double a = 6.0/alpha;
    	double b = r/delta;
    	double b2 = b*b;
    	double b6 = b2*b2*b2;
    	
    	return (epsilon/(1-a))*(a*Math.exp(alpha*(1-b))-(1/b6));
    }
    
    private double du(double r2, double alpha, double epsilon, double delta){
    	
    	double r = Math.sqrt(r2);
    	double a = 6.0/alpha;
    	double b = r/delta;
    	double b2 = b*b;
    	double b6 = b2*b2*b2;
    	
    	return -(6*epsilon/(1-a))*(b*Math.exp(alpha*(1-b))-(1/b6));  	
    }
    
    
	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}
    

    private static final long serialVersionUID = 1L;
    
	protected Boundary boundary;
	protected final double chargeC = ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexCenter];
	protected final double chargeN = ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexN1];
	protected final double chargeP = ConformationNitrogenShellModel.Echarge[SpeciesN2ShellModel.indexP1left];
	protected final double chargeCC, chargeCN, chargeCP;
	protected final double chargeNN, chargeNP;
	protected final double chargePP;
	
	protected final double alpha1 = 10.6; //11.51;
	protected final double alpha2 = 12.7; //10.46;
	protected final double epsilon1 = 36.44; //32.76; // unit K
	protected final double epsilon2 = 35.61; //32.0; // unit K
	protected final double delta1 = 4.531; // unit A
	protected final double delta2 = 3.457; // unit A
	
	protected double rC, r2;
	protected final Vector work, shift;
	protected final Vector[] gradient;

}
