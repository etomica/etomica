/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.IPotentialMolecularTorque;
import etomica.potential.PotentialMolecular;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Tensor3D;
import etomica.units.Kelvin;

/** 
 * P2 6-point potential for Nitrogen.  
 * Reference paper: Etters RD. et al, Phys. Rev B 33(12) 1986
 *
 * 
 * 
 * @author Tai Boon Tan
 */
public class P2Nitrogen extends PotentialMolecular implements IPotentialMolecularTorque {

	public P2Nitrogen(Space space, double rC) {
		super(2, space);
		
		gradient = new Vector[2];
		gradient[0] = space.makeVector();
		gradient[1] = space.makeVector();

		torque = new Vector[2];
	    torque[0] = space.makeVector();
	    torque[1] = space.makeVector();

	    gradientAndTorque = new Vector[2][2];
	    gradientAndTorque[0] = gradient;
	    gradientAndTorque[1] = torque;

	    secDerXr = new Vector[2][3];
	    for(int i=0; i<secDerXr.length; i++){
	    	for(int j=0; j<secDerXr[0].length; j++){
		    	secDerXr[i][j] = space.makeVector();
		    }	
	    }
	    
	    workVec = new Vector[2];
	    for(int i=0; i<workVec.length; i++){
	    	workVec[i] = space.makeVector();
	    }
	    
	    workTorqVec = new Vector[3];
	    for(int i=0; i<workTorqVec.length; i++){
	    	workTorqVec[i] = space.makeVector();
	    }
	    
	    workdTorq = new double[3][3];
	    q = new double[3][3];
	    
	    tensorWork = new DataTensor(space);
	    
		work = space.makeVector();
        duWork = space.makeVector();
        tWork = space.makeVector();
        dr1 = space.makeVector();
        dr2 = space.makeVector();
		shift = space.makeVector();
		com1 = space.makeVector();
		com2 = space.makeVector();
		vectorR = space.makeVector();
		
		C = new double[5];
		// Published values according to the above reference, however found discontinuity 
		//  in the potential
//		C[0] = Kelvin.UNIT.toSim( 415.73107);  //[K]
//		C[1] = Kelvin.UNIT.toSim(-1446.74414); //[KA^-1]
//		C[2] = Kelvin.UNIT.toSim( 2480.73711); //[KA^-2]
//		C[3] = Kelvin.UNIT.toSim(-2766.5419); //[KA^-3]`
//		C[4] = Kelvin.UNIT.toSim( 1574.2809); //[KA^-4]
		
		// Refit the energy to remove discontinuity
//		C[0] = Kelvin.UNIT.toSim(415.7168933551);//415.73107);  //[K]
//		C[1] = Kelvin.UNIT.toSim( -1446.7366463624);//-1446.74414); //[KA^-1]
//		C[2] = Kelvin.UNIT.toSim(2480.7405961429);//2480.73711); //[KA^-2]
//		C[3] = Kelvin.UNIT.toSim(  -2766.5403682242);//-2766.5419); //[KA^-3]`
//		C[4] = Kelvin.UNIT.toSim( 1574.2815730466);//1574.2809); //[KA^-4]
		
		C[0] = Kelvin.UNIT.toSim(  415.7168814243);  //[K]
		C[1] = Kelvin.UNIT.toSim(-1446.7346714270);  //[KA^-1]
		C[2] = Kelvin.UNIT.toSim( 2480.73711);       //[KA^-2]
		C[3] = Kelvin.UNIT.toSim(-2766.5419);        //[KA^-3]`
		C[4] = Kelvin.UNIT.toSim( 1574.2809);        //[KA^-4]
			
        chargeP1P1 = chargeP1 * chargeP1;
        chargeP1P2 = chargeP1 * chargeP2;
        chargeP2P2 = chargeP2 * chargeP2;
        
        this.rC = rC;
	}

	public void setRange(double rC) {
		this.rC = rC;
	}

	public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair){
    	
		double sum = 0.0;
		double r2 = 0.0;

		IMolecule nitrogena = pair.get(0);
		IMolecule nitrogenb = pair.get(1);
		
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
		//System.out.println("<P2Nitrogen> distance: " + Math.sqrt(work.squared()));
		final boolean zeroShift;
		
		if(enablePBC){
			shift.Ea1Tv1(-1,work);
			boundary.nearestImage(work);
			shift.PE(work);
			zeroShift = shift.squared() < 0.1;
		} else {
			zeroShift = true;
		}
		
		r2 = work.squared();
		
		if (r2 > rC*rC){ 
//			System.out.println("TRUNCATED!!!");
//			System.exit(1);
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
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
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
    				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
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
        return sum;																					        
	}
    
	public double virial(IMoleculeList pair) {
		
		IMolecule nitrogena = pair.get(0);
		IMolecule nitrogenb = pair.get(1);
		
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
	    return gradientAndTorque(pair)[0];
	}

    public Vector[][] gradientAndTorque(IMoleculeList pair) {
        IMolecule nitrogena = pair.get(0);
        IMolecule nitrogenb = pair.get(1);

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

        final boolean zeroShift;

        if(enablePBC){
            zeroShift = shift.squared() < 0.1;
        } else {
            zeroShift = true;
        }

        r2 = work.squared();

        gradient[0].E(0.0);
        gradient[1].E(0.0);

        torque[0].E(0);
        torque[1].E(0);

        if (r2 > rC*rC){ 
            //System.out.println("TRUNCATED!!!");
            return gradientAndTorque;
        }
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
             *  non-bonded atoms between the 2 molecules 
             */
            
            for (int i=0; i<2; i++){
                dr1.Ev1Mv2(nitrogena.getChildList().get(i).getPosition(), com1);
                for (int j=0; j<2; j++){
                    work.Ev1Mv2(nitrogena.getChildList().get(i).getPosition(), nitrogenb.getChildList().get(j).getPosition());
                    dr2.Ev1Mv2(nitrogenb.getChildList().get(j).getPosition(), com2);
                    r2 = work.squared();
                    
                    if(Math.sqrt(r2) >= R1){            // R >= R1
                        duWork.Ea1Tv1(dURgtR1(r2)/r2, work);
                    } else if (Math.sqrt(r2) < R1 && Math.sqrt(r2) >= R0){  // R1 > R >= R0
                        duWork.Ea1Tv1(dUR1gtRgteqR0(r2)/r2, work);
                    } else if (Math.sqrt(r2) < R0){     // R < R0
                        duWork.Ea1Tv1(dURltR0(r2)/r2, work);
                    }
                    gradient[0].PE(duWork);
                    tWork.E(dr1);
                    tWork.XE(duWork);
                    torque[0].ME(tWork);
                    tWork.E(dr2);
                    tWork.XE(duWork);
                    torque[1].PE(tWork);

                }
            }
                        
            doGradientNoShift(Pa1l, Pb1l, chargeP1P1);
            doGradientNoShift(Pa2l, Pb1l, chargeP1P2);
            doGradientNoShift(Pa2r, Pb1l, chargeP1P2);
            doGradientNoShift(Pa1r, Pb1l, chargeP1P1);
            
            doGradientNoShift(Pa1l, Pb2l, chargeP1P2);
            doGradientNoShift(Pa2l, Pb2l, chargeP2P2);
            doGradientNoShift(Pa2r, Pb2l, chargeP2P2);
            doGradientNoShift(Pa1r, Pb2l, chargeP1P2);
            
            doGradientNoShift(Pa1l, Pb2r, chargeP1P2);
            doGradientNoShift(Pa2l, Pb2r, chargeP2P2);
            doGradientNoShift(Pa2r, Pb2r, chargeP2P2);
            doGradientNoShift(Pa1r, Pb2r, chargeP1P2);

   		 
			
            doGradientNoShift(Pa1l, Pb1r, chargeP1P1);
            doGradientNoShift(Pa2l, Pb1r, chargeP1P2);
            doGradientNoShift(Pa2r, Pb1r, chargeP1P2);
            doGradientNoShift(Pa1r, Pb1r, chargeP1P1);
        }
        else {
            
            /*
             * 'for' loop for 4 pairs van der Waals interaction between the 
             *  non-bonded atoms between the 2 molecules 
             * 
             * 
             */
            
            for (int i=0; i<2; i++){
                dr2.Ev1Mv2(nitrogenb.getChildList().get(i).getPosition(), com2);
                Vector dist = (nitrogenb.getChildList().get(i)).getPosition();
                shift.TE(-1.0);
                shift.PE(dist);
                
                for (int j=0; j<2; j++){
                    dr1.Ev1Mv2(nitrogena.getChildList().get(j).getPosition(), com1);
                    work.Ev1Mv2(nitrogena.getChildList().get(j).getPosition(), shift);
                    r2 = work.squared();
                    
                    if(Math.sqrt(r2) >= R1){            // R >= R1
                        duWork.Ea1Tv1(dURgtR1(r2)/r2, work);
                    } else if (Math.sqrt(r2) < R1 && Math.sqrt(r2) >= R0){  // R1 > R >= R0
                        duWork.Ea1Tv1(dUR1gtRgteqR0(r2)/r2, work);
                    } else if (Math.sqrt(r2) < R0){     // R < R0
                        duWork.Ea1Tv1(dURltR0(r2)/r2, work);
                    }
                    
                    gradient[0].PE(duWork);
                    tWork.E(dr1);
                    tWork.XE(duWork);
                    torque[0].ME(tWork);
                    tWork.E(dr2);
                    tWork.XE(duWork);
                    torque[1].PE(tWork);
                }
                shift.ME(dist);
                shift.TE(-1.0);
            }
            
            doGradientShift(Pa1l, Pb1l, chargeP1P1);
            doGradientShift(Pa2l, Pb1l, chargeP1P2);
            doGradientShift(Pa2r, Pb1l, chargeP1P2);
            doGradientShift(Pa1r, Pb1l, chargeP1P1);
            
            doGradientShift(Pa1l, Pb2l, chargeP1P2);
            doGradientShift(Pa2l, Pb2l, chargeP2P2);
            doGradientShift(Pa2r, Pb2l, chargeP2P2);
            doGradientShift(Pa1r, Pb2l, chargeP1P2);
            
            doGradientShift(Pa1l, Pb2r, chargeP1P2);
            doGradientShift(Pa2l, Pb2r, chargeP2P2);
            doGradientShift(Pa2r, Pb2r, chargeP2P2);
            doGradientShift(Pa1r, Pb2r, chargeP1P2);

            doGradientShift(Pa1l, Pb1r, chargeP1P1);
            doGradientShift(Pa2l, Pb1r, chargeP1P2);
            doGradientShift(Pa2r, Pb1r, chargeP1P2);
            doGradientShift(Pa1r, Pb1r, chargeP1P1);
        }
        
        gradient[1].Ea1Tv1(-1, gradient[0]);
      
        return gradientAndTorque;
    }

    protected void doGradientShift(Vector v1, Vector v2, double cc) {
        shift.TE(-1.0);
        shift.PE(v2);
        work.Ev1Mv2(v1, shift);
        r2 = work.squared();
        shift.ME(v2);
        duWork.Ea1Tv1(-cc/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(duWork);
        dr1.Ev1Mv2(v1, com1);
        tWork.E(dr1);
        tWork.XE(duWork);
        torque[0].ME(tWork);
        dr2.Ev1Mv2(v2, com2);
        tWork.E(dr2);
        tWork.XE(duWork);
        torque[1].PE(tWork);
        shift.TE(-1.0);
    }

    protected void doGradientNoShift(Vector v1, Vector v2, double cc) {
        work.Ev1Mv2(v1, v2);
        r2 = work.squared();
        duWork.Ea1Tv1(-cc/(r2*Math.sqrt(r2)), work);
        gradient[0].PE(duWork);
        dr1.Ev1Mv2(v2, com1);
        tWork.E(dr1);
        tWork.XE(duWork);
        torque[0].ME(tWork);
        dr2.Ev1Mv2(v2, com2);
        tWork.E(dr2);
        tWork.XE(duWork);
        torque[1].PE(tWork);
    }

    public DataTensor secondDerivative(IMoleculeList pair){
    	
    	DataTensor tensor = new DataTensor(space);
    	DataTensor sumTensor = new DataTensor(space);
    	
    	double r;
		double r2 = 0.0;

		IMolecule nitrogena = pair.get(0);
		IMolecule nitrogenb = pair.get(1);
		
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
		//System.out.println("<P2Nitrogen> distance: " + Math.sqrt(work.squared()));
		final boolean zeroShift;
		
		if(enablePBC){
			shift.Ea1Tv1(-1,work);
			boundary.nearestImage(work);
			shift.PE(work);
			zeroShift = shift.squared() < 0.1;
		} else {
			zeroShift = true;
		}
		
		r2 = work.squared();
		
		if (r2 > rC*rC){ 
			//System.out.println("TRUNCATED!!!");
			return sumTensor;
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
    				
    				vectorR.Ev1Mv2(dist, (nitrogenb.getChildList().get(j)).getPosition());
    				double distr2 = vectorR.squared();
    				tensor.x.Ev1v2(vectorR, vectorR);
    				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
    					tensor.TE(1.0/(distr2*distr2)*(dURgtR1(distr2) - d2URgtR1(distr2)));
    					tensor.x.PEa1Tt1(-dURgtR1(distr2)/distr2, identity);
    					sumTensor.PE(tensor);
    				
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					tensor.TE(1.0/(distr2*distr2)*(dUR1gtRgteqR0(distr2) - d2UR1gtRgteqR0(distr2)));
    					tensor.x.PEa1Tt1(-dUR1gtRgteqR0(distr2)/distr2, identity);
    					sumTensor.PE(tensor);
    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					tensor.TE(1.0/(distr2*distr2)*(dURltR0(distr2) - d2URltR0(distr2)));
    					tensor.x.PEa1Tt1(-dURltR0(distr2)/distr2, identity);
    					sumTensor.PE(tensor);
    				}

    			}

    		}
    		        	
        	//Pa1l
    		vectorR.Ev1Mv2(Pa1l, Pb1l);
			r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
			
			
            vectorR.Ev1Mv2(Pa1l, Pb2l);
        	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
			
			
            vectorR.Ev1Mv2(Pa1l, Pb2r);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
            vectorR.Ev1Mv2(Pa1l, Pb1r);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
        
            //Pa2l
            vectorR.Ev1Mv2(Pa2l, Pb1l);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);

            vectorR.Ev1Mv2(Pa2l, Pb2l);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
			vectorR.Ev1Mv2(Pa2l, Pb2r);
		   	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
			
            
			vectorR.Ev1Mv2(Pa2l, Pb1r);
		   	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
		
            //Pa2r
            vectorR.Ev1Mv2(Pa2r, Pb1l);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
            vectorR.Ev1Mv2(Pa2r, Pb2l);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
            vectorR.Ev1Mv2(Pa2r, Pb2r);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
            vectorR.Ev1Mv2(Pa2r, Pb1r);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
       
            //Pa1r
            vectorR.Ev1Mv2(Pa1r, Pb1l);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
            
            vectorR.Ev1Mv2(Pa1r, Pb2l);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
            
            vectorR.Ev1Mv2(Pa1r, Pb2r);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
            
            vectorR.Ev1Mv2(Pa1r, Pb1r);
           	r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
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
    				
    				vectorR.Ev1Mv2((nitrogena.getChildList().get(j)).getPosition(), shift);
    				double distr2 = vectorR.squared();
    				tensor.x.Ev1v2(vectorR, vectorR);
    				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
    					tensor.TE(1.0/(distr2*distr2)*(dURgtR1(distr2) - d2URgtR1(distr2)));
    					tensor.x.PEa1Tt1(-dURgtR1(distr2)/distr2, identity);
    					sumTensor.PE(tensor);
    				
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					tensor.TE(1.0/(distr2*distr2)*(dUR1gtRgteqR0(distr2) - d2UR1gtRgteqR0(distr2)));
    					tensor.x.PEa1Tt1(-dUR1gtRgteqR0(distr2)/distr2, identity);
    					sumTensor.PE(tensor);
    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					tensor.TE(1.0/(distr2*distr2)*(dURltR0(distr2) - d2URltR0(distr2)));
    					tensor.x.PEa1Tt1(-dURltR0(distr2)/distr2, identity);
    					sumTensor.PE(tensor);
    				}
    				    				
    			}
    			shift.ME(dist);
    			shift.TE(-1.0);
    		}
    		
        	shift.TE(-1.0);
        	shift.PE(Pb1l);
            vectorR.Ev1Mv2(Pa1l, shift);
            shift.ME(Pb1l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);      
            vectorR.Ev1Mv2(Pa2l, shift);
            shift.ME(Pb1l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);
            vectorR.Ev1Mv2(Pa2r, shift);
            shift.ME(Pb1l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1l);
            vectorR.Ev1Mv2(Pa1r, shift);
            shift.ME(Pb1l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);
            
            ////////////
        	shift.TE(-1.0);
            shift.PE(Pb2l);
            vectorR.Ev1Mv2(Pa1l, shift);
            shift.ME(Pb2l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            vectorR.Ev1Mv2(Pa2l, shift);
            shift.ME(Pb2l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            vectorR.Ev1Mv2(Pa2r, shift);
            shift.ME(Pb2l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);
            

        	shift.TE(-1.0);
            shift.PE(Pb2l);
            vectorR.Ev1Mv2(Pa1r, shift);
            shift.ME(Pb2l);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);
            
            //////////////////////
        	shift.TE(-1.0);
            shift.PE(Pb2r);
            vectorR.Ev1Mv2(Pa1l, shift);
            shift.ME(Pb2r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            vectorR.Ev1Mv2(Pa2l, shift);
            shift.ME(Pb2r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            vectorR.Ev1Mv2(Pa2r, shift);
            shift.ME(Pb2r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP2P2/r) - (2*chargeP2P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP2P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb2r);
            vectorR.Ev1Mv2(Pa1r, shift);
            shift.ME(Pb2r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);
            
            /////////////
        	shift.TE(-1.0);
            shift.PE(Pb1r);
            vectorR.Ev1Mv2(Pa1l, shift);
            shift.ME(Pb1r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            vectorR.Ev1Mv2(Pa2l, shift);
            shift.ME(Pb1r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            vectorR.Ev1Mv2(Pa2r, shift);
            shift.ME(Pb1r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P2/r) - (2*chargeP1P2/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P2/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);

        	shift.TE(-1.0);
            shift.PE(Pb1r);
            vectorR.Ev1Mv2(Pa1r, shift);
            shift.ME(Pb1r);
            r2 = vectorR.squared();
			r = Math.sqrt(r2);
			tensor.x.Ev1v2(vectorR, vectorR);
			tensor.TE(1.0/(r2*r2)*( (-chargeP1P1/r) - (2*chargeP1P1/r)));
			tensor.x.PEa1Tt1(-(-chargeP1P1/r)/r2, identity);
			sumTensor.PE(tensor);
        	shift.TE(-1.0);
            
        }
        return sumTensor;																					        
	}
	
    
    public Vector[][] secondDerivativeXr(IMoleculeList pair){
    	
    	DataTensor tensor = new DataTensor(space);
    	
    	//  secDerXr[0] returns quantity for (molA - molB)
    	//  secDerXr[1] returns quantity for (molB - molA)
    	for(int i=0; i<secDerXr.length; i++){
  	    	for(int j=0; j<secDerXr[0].length; j++){
  		    	secDerXr[i][j].E(0.0);
  		    }	
  	    }
    	
		double r2 = 0.0;

		IMolecule nitrogena = pair.get(0);
		IMolecule nitrogenb = pair.get(1);
		
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

		final boolean zeroShift;
		
		if(enablePBC){
			zeroShift = shift.squared() < 0.1;
		} else {
			zeroShift = true;
		}
		
		r2 = work.squared();
		
		if (r2 > rC*rC){ 
			return secDerXr;
		}
		
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
    			dr1.Ev1Mv2(nitrogena.getChildList().get(i).getPosition(), com1);
    			
    			for (int j=0; j<2; j++){
    				
    				vectorR.Ev1Mv2((nitrogena.getChildList().get(i)).getPosition(), (nitrogenb.getChildList().get(j)).getPosition());
    				dr2.Ev1Mv2(nitrogenb.getChildList().get(j).getPosition(), com2);
    				
    				double distr2 = vectorR.squared();
    				tensor.x.Ev1v2(vectorR, vectorR);
    				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
    					tensor.TE(1.0/(distr2*distr2)*(dURgtR1(distr2) - d2URgtR1(distr2)));
    					tensor.x.PEa1Tt1(-dURgtR1(distr2)/distr2, identity);
    				    					
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					tensor.TE(1.0/(distr2*distr2)*(dUR1gtRgteqR0(distr2) - d2UR1gtRgteqR0(distr2)));
    					tensor.x.PEa1Tt1(-dUR1gtRgteqR0(distr2)/distr2, identity);
    					    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					tensor.TE(1.0/(distr2*distr2)*(dURltR0(distr2) - d2URltR0(distr2)));
    					tensor.x.PEa1Tt1(-dURltR0(distr2)/distr2, identity);
    					
    				}
    				
    				for(int iloop=0; iloop<3; iloop++){
						workVec[0].E(new double[]{tensor.x.component(iloop, 0),
												  tensor.x.component(iloop, 1), 
								                  tensor.x.component(iloop, 2)});
						workVec[1].E(workVec[0]);
						workVec[0].XE(dr1);
						workVec[1].XE(dr2);
						
						secDerXr[0][iloop].PE(workVec[0]);
						secDerXr[1][iloop].PE(workVec[1]);
					}
    			}
    		}
    		        	
        	//Pa1l
    		doSecDerXrNoShift(Pa1l, Pb1l, chargeP1P1);
    		doSecDerXrNoShift(Pa1l, Pb2l, chargeP1P2);
        	doSecDerXrNoShift(Pa1l, Pb2r, chargeP1P2);
           	doSecDerXrNoShift(Pa1l, Pb1r, chargeP1P1);
        
            //Pa2l
           	doSecDerXrNoShift(Pa2l, Pb1l, chargeP1P2);
           	doSecDerXrNoShift(Pa2l, Pb2l, chargeP2P2);
           	doSecDerXrNoShift(Pa2l, Pb2r, chargeP2P2);
		   	doSecDerXrNoShift(Pa2l, Pb1r, chargeP1P2);
		
            //Pa2r
		   	doSecDerXrNoShift(Pa2r, Pb1l, chargeP1P2);
           	doSecDerXrNoShift(Pa2r, Pb2l, chargeP2P2);
           	doSecDerXrNoShift(Pa2r, Pb2r, chargeP2P2);
           	doSecDerXrNoShift(Pa2r, Pb1r, chargeP1P2);
        
            //Pa1r
           	doSecDerXrNoShift(Pa1r, Pb1l, chargeP1P1);
           	doSecDerXrNoShift(Pa1r, Pb2l, chargeP1P2);
           	doSecDerXrNoShift(Pa1r, Pb2r, chargeP1P2);
           	doSecDerXrNoShift(Pa1r, Pb1r, chargeP1P1);
           	
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
    			
    			dr2.Ev1Mv2(nitrogenb.getChildList().get(i).getPosition(), com2);
    			
    			for (int j=0; j<2; j++){
    				
    				dr1.Ev1Mv2(nitrogena.getChildList().get(j).getPosition(), com1);
    				
    				vectorR.Ev1Mv2((nitrogena.getChildList().get(j)).getPosition(), shift);
    				double distr2 = vectorR.squared();
    				tensor.x.Ev1v2(vectorR, vectorR);
    				
        				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
    					tensor.TE(1.0/(distr2*distr2)*(dURgtR1(distr2) - d2URgtR1(distr2)));
    					tensor.x.PEa1Tt1(-dURgtR1(distr2)/distr2, identity);
    				
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					tensor.TE(1.0/(distr2*distr2)*(dUR1gtRgteqR0(distr2) - d2UR1gtRgteqR0(distr2)));
    					tensor.x.PEa1Tt1(-dUR1gtRgteqR0(distr2)/distr2, identity);
    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					tensor.TE(1.0/(distr2*distr2)*(dURltR0(distr2) - d2URltR0(distr2)));
    					tensor.x.PEa1Tt1(-dURltR0(distr2)/distr2, identity);
    					
    				}
    				
    				for(int iloop=0; iloop<3; iloop++){
						workVec[0].E(new double[]{tensor.x.component(iloop, 0),
												  tensor.x.component(iloop, 1), 
								                  tensor.x.component(iloop, 2)});
						workVec[1].E(workVec[0]);
						workVec[0].XE(dr1);
						workVec[1].XE(dr2);
						
						secDerXr[0][iloop].PE(workVec[0]);
						secDerXr[1][iloop].PE(workVec[1]);
					}				
    			}
    			shift.ME(dist);
    			shift.TE(-1.0);
    		}
    		
        	doSecDerXrWShift(Pa1l, Pb1l, chargeP1P1);
            doSecDerXrWShift(Pa2l, Pb1l, chargeP1P2);
            doSecDerXrWShift(Pa2r, Pb1l, chargeP1P2);
            doSecDerXrWShift(Pa1r, Pb1l, chargeP1P1);
            
            ////////////
            doSecDerXrWShift(Pa1l, Pb2l, chargeP1P2);
            doSecDerXrWShift(Pa2l, Pb2l, chargeP2P2);
            doSecDerXrWShift(Pa2r, Pb2l, chargeP2P2);
            doSecDerXrWShift(Pa1r, Pb2l, chargeP1P2);
            
            //////////////////////
        	doSecDerXrWShift(Pa1l, Pb2r, chargeP1P2);
            doSecDerXrWShift(Pa2l, Pb2r, chargeP2P2);
            doSecDerXrWShift(Pa2r, Pb2r, chargeP2P2);
            doSecDerXrWShift(Pa1r, Pb2r, chargeP1P2);
            
            /////////////
        	doSecDerXrWShift(Pa1l, Pb1r, chargeP1P1);
            doSecDerXrWShift(Pa2l, Pb1r, chargeP1P2);
            doSecDerXrWShift(Pa2r, Pb1r, chargeP1P2);
            doSecDerXrWShift(Pa1r, Pb1r, chargeP1P1);
            
        }
        
        return secDerXr;																					        
	}
    
    public void doSecDerXrNoShift(Vector v1, Vector v2, double cc){
    	vectorR.Ev1Mv2(v1, v2);
		r2 = vectorR.squared();
		double r = Math.sqrt(r2);
		tensorWork.x.Ev1v2(vectorR, vectorR);
		tensorWork.TE(1.0/(r2*r2)*( (-cc/r) - (2*cc/r)));
		tensorWork.x.PEa1Tt1(-(-cc/r)/r2, identity);
		
		dr1.Ev1Mv2(v1, com1);
		dr2.Ev1Mv2(v2, com2);
		
		for(int iloop=0; iloop<3; iloop++){
			workVec[0].E(new double[]{tensorWork.x.component(iloop, 0),
									  tensorWork.x.component(iloop, 1), 
					                  tensorWork.x.component(iloop, 2)});
			workVec[1].E(workVec[0]);
			workVec[0].XE(dr1);
			workVec[1].XE(dr2);
			
			secDerXr[0][iloop].PE(workVec[0]);
			secDerXr[1][iloop].PE(workVec[1]);
		}
    	
    }
    
    public void doSecDerXrWShift(Vector v1, Vector v2, double cc){
    	shift.TE(-1.0);
    	shift.PE(v2);
        vectorR.Ev1Mv2(v1, shift);
        shift.ME(v2);
        shift.TE(-1.0);
        
        r2 = vectorR.squared();
		double r = Math.sqrt(r2);
		tensorWork.x.Ev1v2(vectorR, vectorR);
		tensorWork.TE(1.0/(r2*r2)*( (-cc/r) - (2*cc/r)));
		tensorWork.x.PEa1Tt1(-(-cc/r)/r2, identity);
		
		dr1.Ev1Mv2(v1, com1);
		dr2.Ev1Mv2(v2, com2);
		
		for(int iloop=0; iloop<3; iloop++){
			workVec[0].E(new double[]{tensorWork.x.component(iloop, 0),
									  tensorWork.x.component(iloop, 1), 
					                  tensorWork.x.component(iloop, 2)});
			workVec[1].E(workVec[0]);
			workVec[0].XE(dr1);
			workVec[1].XE(dr2);
			
			secDerXr[0][iloop].PE(workVec[0]);
			secDerXr[1][iloop].PE(workVec[1]);
		}
    	
    }
    
    
    public Tensor secondDerivativeXrRotRot(IMoleculeList pair){
    	
    	DataTensor tensor = new DataTensor(space);
    	Tensor tensorq = space.makeTensor(); 
    	
    	//  secDerXr[0] returns quantity for (molA - molB)
    	//  secDerXr[1] returns quantity for (molB - molA)
    	for(int i=0; i<secDerXr.length; i++){
  	    	for(int j=0; j<secDerXr[0].length; j++){
  		    	secDerXr[i][j].E(0.0);
  		    }	
  	    }
    	
    	for(int i=0; i<q.length; i++){
  	    	for(int j=0; j<q[0].length; j++){
  		    	q[i][j] = 0.0;
  		    }	
  	    }
    	
		double r2 = 0.0;

		IMolecule nitrogena = pair.get(0);
		IMolecule nitrogenb = pair.get(1);
		
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

		final boolean zeroShift;
		
		if(enablePBC){
			zeroShift = shift.squared() < 0.1;
		} else {
			zeroShift = true;
		}
		
		r2 = work.squared();
		
		if (r2 > rC*rC){ 
			return tensorq;
		}
		
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
    			dr1.Ev1Mv2(nitrogena.getChildList().get(i).getPosition(), com1);
    			
    			for(int iTorq=0; iTorq<3; iTorq++){ //iTorq = row
    	  	    	for(int jTorq=0; jTorq<3; jTorq++){ //jTorq = column
    	  		    	workdTorq[iTorq][jTorq] = 0.0;
    	  		    }	
    	  	    }
    			
    			for (int j=0; j<2; j++){
    				
    				vectorR.Ev1Mv2((nitrogena.getChildList().get(i)).getPosition(), (nitrogenb.getChildList().get(j)).getPosition());
    				dr2.Ev1Mv2(nitrogenb.getChildList().get(j).getPosition(), com2);
    				
    				double distr2 = vectorR.squared();
    				tensor.x.Ev1v2(vectorR, vectorR);
    				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
    					tensor.TE(1.0/(distr2*distr2)*(dURgtR1(distr2) - d2URgtR1(distr2)));
    					tensor.x.PEa1Tt1(-dURgtR1(distr2)/distr2, identity);
    				    					
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					tensor.TE(1.0/(distr2*distr2)*(dUR1gtRgteqR0(distr2) - d2UR1gtRgteqR0(distr2)));
    					tensor.x.PEa1Tt1(-dUR1gtRgteqR0(distr2)/distr2, identity);
    					    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					tensor.TE(1.0/(distr2*distr2)*(dURltR0(distr2) - d2URltR0(distr2)));
    					tensor.x.PEa1Tt1(-dURltR0(distr2)/distr2, identity);
    					
    				}

    				for(int iloop=0; iloop<3; iloop++){
						workVec[0].E(new double[]{tensor.x.component(iloop, 0),
												  tensor.x.component(iloop, 1), 
								                  tensor.x.component(iloop, 2)});
						workVec[1].E(workVec[0]);
						//workVec[0].XE(dr1);
						workVec[1].XE(dr2);
						
						workdTorq[0][iloop] += workVec[1].getX(0);
						workdTorq[1][iloop] += workVec[1].getX(1);
						workdTorq[2][iloop] += workVec[1].getX(2);
						
					}
    				
    			}
    			
    			for(int iq=0; iq<3; iq++){
    				
    				tWork.E(workdTorq[iq]);
	    			tWork.XE(dr1);
	    			for(int jq=0; jq<3; jq++){
	    				q[jq][iq] += tWork.getX(jq); 
	    			}
    			}
    			
    		}
    		        	
        	//Pa1l
    		doSecDerXrRotRotNoShift(Pa1l, Pb1l, chargeP1P1);
    		doSecDerXrRotRotNoShift(Pa1l, Pb2l, chargeP1P2);
        	doSecDerXrRotRotNoShift(Pa1l, Pb2r, chargeP1P2);
           	doSecDerXrRotRotNoShift(Pa1l, Pb1r, chargeP1P1);
        
            //Pa2l
           	doSecDerXrRotRotNoShift(Pa2l, Pb1l, chargeP1P2);
           	doSecDerXrRotRotNoShift(Pa2l, Pb2l, chargeP2P2);
           	doSecDerXrRotRotNoShift(Pa2l, Pb2r, chargeP2P2);
		   	doSecDerXrRotRotNoShift(Pa2l, Pb1r, chargeP1P2);
		
            //Pa2r
		   	doSecDerXrRotRotNoShift(Pa2r, Pb1l, chargeP1P2);
           	doSecDerXrRotRotNoShift(Pa2r, Pb2l, chargeP2P2);
           	doSecDerXrRotRotNoShift(Pa2r, Pb2r, chargeP2P2);
           	doSecDerXrRotRotNoShift(Pa2r, Pb1r, chargeP1P2);
        
            //Pa1r
           	doSecDerXrRotRotNoShift(Pa1r, Pb1l, chargeP1P1);
           	doSecDerXrRotRotNoShift(Pa1r, Pb2l, chargeP1P2);
           	doSecDerXrRotRotNoShift(Pa1r, Pb2r, chargeP1P2);
           	doSecDerXrRotRotNoShift(Pa1r, Pb1r, chargeP1P1);
           	
        } 
        
        else {
        	
    		/*
    		 * 'for' loop for 4 pairs van der Waals interaction between the 
    		 * 	non-bonded atoms between the 2 molecules 
    		 * 
    		 * 
    		 */
        	
    		for (int i=0; i<2; i++){
    			Vector dist = (nitrogena.getChildList().get(i)).getPosition();
    			shift.PE(dist);
    			
    			dr1.Ev1Mv2(nitrogena.getChildList().get(i).getPosition(), com1);
    			
    			for(int iTorq=0; iTorq<3; iTorq++){ //iTorq = row
    	  	    	for(int jTorq=0; jTorq<3; jTorq++){ //jTorq = column
    	  		    	workdTorq[iTorq][jTorq] = 0.0;
    	  		    }	
    	  	    }
    			
    			for (int j=0; j<2; j++){
    				
    				vectorR.Ev1Mv2((nitrogenb.getChildList().get(j)).getPosition(), shift);
    				double distr2 = vectorR.squared();
    				tensor.x.Ev1v2(vectorR, vectorR);
    				dr2.Ev1Mv2(nitrogenb.getChildList().get(j).getPosition(), com2);
        				
    				if(Math.sqrt(distr2) >= R1){            // R >= R1
    					tensor.TE(1.0/(distr2*distr2)*(dURgtR1(distr2) - d2URgtR1(distr2)));
    					tensor.x.PEa1Tt1(-dURgtR1(distr2)/distr2, identity);
    				
    				} else if (Math.sqrt(distr2) < R1 && Math.sqrt(distr2) >= R0){  // R1 > R >= R0
    					tensor.TE(1.0/(distr2*distr2)*(dUR1gtRgteqR0(distr2) - d2UR1gtRgteqR0(distr2)));
    					tensor.x.PEa1Tt1(-dUR1gtRgteqR0(distr2)/distr2, identity);
    					
    				} else if (Math.sqrt(distr2) < R0){   	// R < R0
    					tensor.TE(1.0/(distr2*distr2)*(dURltR0(distr2) - d2URltR0(distr2)));
    					tensor.x.PEa1Tt1(-dURltR0(distr2)/distr2, identity);
    					
    				}
    				
    				for(int iloop=0; iloop<3; iloop++){
						workVec[0].E(new double[]{tensor.x.component(iloop, 0),
												  tensor.x.component(iloop, 1), 
								                  tensor.x.component(iloop, 2)});
						workVec[1].E(workVec[0]);
						//workVec[0].XE(dr1);
						workVec[1].XE(dr2);
						
						workdTorq[0][iloop] += workVec[1].getX(0);
						workdTorq[1][iloop] += workVec[1].getX(1);
						workdTorq[2][iloop] += workVec[1].getX(2);
					}		
    			}
    			
    			for(int iq=0; iq<3; iq++){	
    				tWork.E(workdTorq[iq]);
	    			tWork.XE(dr1);
	    			for(int jq=0; jq<3; jq++){
	    				q[jq][iq] += tWork.getX(jq); 
	    			}
    			}
    			
    			shift.ME(dist);
    		}
    		
        	doSecDerXrRotRotWShift(Pa1l, Pb1l, chargeP1P1);
            doSecDerXrRotRotWShift(Pa2l, Pb1l, chargeP1P2);
            doSecDerXrRotRotWShift(Pa2r, Pb1l, chargeP1P2);
            doSecDerXrRotRotWShift(Pa1r, Pb1l, chargeP1P1);
            
            ////////////
            doSecDerXrRotRotWShift(Pa1l, Pb2l, chargeP1P2);
            doSecDerXrRotRotWShift(Pa2l, Pb2l, chargeP2P2);
            doSecDerXrRotRotWShift(Pa2r, Pb2l, chargeP2P2);
            doSecDerXrRotRotWShift(Pa1r, Pb2l, chargeP1P2);
            
            //////////////////////
        	doSecDerXrRotRotWShift(Pa1l, Pb2r, chargeP1P2);
            doSecDerXrRotRotWShift(Pa2l, Pb2r, chargeP2P2);
            doSecDerXrRotRotWShift(Pa2r, Pb2r, chargeP2P2);
            doSecDerXrRotRotWShift(Pa1r, Pb2r, chargeP1P2);
            
            /////////////
        	doSecDerXrRotRotWShift(Pa1l, Pb1r, chargeP1P1);
            doSecDerXrRotRotWShift(Pa2l, Pb1r, chargeP1P2);
            doSecDerXrRotRotWShift(Pa2r, Pb1r, chargeP1P2);
            doSecDerXrRotRotWShift(Pa1r, Pb1r, chargeP1P1);
            
        }
        
        tensorq.E(q);
        
        return tensorq;																					        
	
    }
    
    public void doSecDerXrRotRotNoShift(Vector v1, Vector v2, double cc){
    	vectorR.Ev1Mv2(v1, v2);
		r2 = vectorR.squared();
		double r = Math.sqrt(r2);
		tensorWork.x.Ev1v2(vectorR, vectorR);
		tensorWork.TE(1.0/(r2*r2)*( (-cc/r) - (2*cc/r)));
		tensorWork.x.PEa1Tt1(-(-cc/r)/r2, identity);
		
		dr1.Ev1Mv2(v1, com1);
		dr2.Ev1Mv2(v2, com2);
		
		for(int iloop=0; iloop<3; iloop++){
			workVec[0].E(new double[]{tensorWork.x.component(iloop, 0),
									  tensorWork.x.component(iloop, 1), 
					                  tensorWork.x.component(iloop, 2)});
			workVec[1].E(workVec[0]);
			//workVec[0].XE(dr1);
			workVec[1].XE(dr2);
			
			workdTorq[0][iloop] = workVec[1].getX(0);
			workdTorq[1][iloop] = workVec[1].getX(1);
			workdTorq[2][iloop] = workVec[1].getX(2);			
		}
		
		for(int iq=0; iq<3; iq++){	
			tWork.E(workdTorq[iq]);
			tWork.XE(dr1);
			for(int jq=0; jq<3; jq++){
				q[jq][iq] += tWork.getX(jq); 
			}
		}
    	
    }
    
    public void doSecDerXrRotRotWShift(Vector v1, Vector v2, double cc){
    	shift.TE(-1.0);
    	shift.PE(v2);
        vectorR.Ev1Mv2(v1, shift);
        shift.ME(v2);
        shift.TE(-1.0);
        
        r2 = vectorR.squared();
		double r = Math.sqrt(r2);
		tensorWork.x.Ev1v2(vectorR, vectorR);
		tensorWork.TE(1.0/(r2*r2)*( (-cc/r) - (2*cc/r)));
		tensorWork.x.PEa1Tt1(-(-cc/r)/r2, identity);
		
		dr1.Ev1Mv2(v1, com1);
		dr2.Ev1Mv2(v2, com2);
		
		for(int iloop=0; iloop<3; iloop++){
			workVec[0].E(new double[]{tensorWork.x.component(iloop, 0),
									  tensorWork.x.component(iloop, 1), 
					                  tensorWork.x.component(iloop, 2)});
			workVec[1].E(workVec[0]);
			//workVec[0].XE(dr1);
			workVec[1].XE(dr2);
			
			workdTorq[0][iloop] = workVec[1].getX(0);
			workdTorq[1][iloop] = workVec[1].getX(1);
			workdTorq[2][iloop] = workVec[1].getX(2);
		}
		
		for(int iq=0; iq<3; iq++){	
			tWork.E(workdTorq[iq]);
			tWork.XE(dr1);
			for(int jq=0; jq<3; jq++){
				q[jq][iq] += tWork.getX(jq); 
			}
		}
    	
    }
    
	public Vector[] gradient(IMoleculeList atoms, Tensor pressureTensor) {
		// TODO Auto-generated method stub
		return null;
	}
	/*
	 * energy potential
	 */
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
    /*
     * first derivative
     */
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

    /*
     * second derivative
     */
    
    private double d2URgtR1(double r2){
    	double r = Math.sqrt(r2);
    	return  (alpha1*alpha1*r*r)*A1*Math.exp(-alpha1*r) - (6*7)*B1/(r2*r2*r2);
    
    }
    
    private double d2UR1gtRgteqR0(double r2){
    	double r = Math.sqrt(r2);
    	double rDiff = (r-R0); 
    	double sumU = 12*C[4]*rDiff*rDiff + 6*C[3]*rDiff + 2*C[2];
    	
    	return r*r*sumU - (6*7)*B1/(r2*r2*r2);
    }
    
    private double d2URltR0(double r2){
    	double r = Math.sqrt(r2);
    	return (alpha2*alpha2*r*r)*A2*Math.exp(-alpha2*r) - (6*7)*B1/(r2*r2*r2);
    	
    }
    
	public boolean isEnablePBC() {
		return enablePBC;
	}

	public void setEnablePBC(boolean enablePBC) {
		this.enablePBC = enablePBC;
	}
    
    
    public double getRange() {
        return rC;
    }
    
    private static final long serialVersionUID = 1L;
        
    protected final Vector[] gradient;
    protected final Vector[] torque;
    protected final Vector[] workVec;
    protected final Vector[] workTorqVec;
    protected final double [][] workdTorq, q;
    protected final Vector[][] gradientAndTorque, secDerXr;
	protected Boundary boundary;
	protected DataTensor tensorWork;
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
	
	protected final Vector work, shift, duWork, tWork, dr1, dr2;
	protected final Vector com1, com2, vectorR;
	protected double rC, r2;
	protected boolean enablePBC = true;
	protected final Tensor identity = new Tensor3D(new double[][]{{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}});
}
