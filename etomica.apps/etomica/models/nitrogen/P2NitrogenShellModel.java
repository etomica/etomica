
package etomica.models.nitrogen;

import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.potential.PotentialMolecular;
import etomica.space.ISpace;

/** 
 * P2 shell-model potential for Nitrogen.  
 *  Reference: 1. Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
 *              	phases, JCP 112(15) 6745 (2000)
 * 			   2. Jordan P.C. et al, Towards phase transferable potential functions: Methodology
 * 					and application to nitrogen, JCP 103(6) 2272 (1995)
 * @author Tai Boon Tan
 */
public class P2NitrogenShellModel extends PotentialMolecular {

	public P2NitrogenShellModel(ISpace space) {
		super(2, space);
		work = space.makeVector();
		shift = space.makeVector();
		
        chargeCC = chargeC * chargeC;
        chargeCN = chargeC * chargeN;
        chargeCP = chargeC * chargeP;
        
        chargeNN = chargeN * chargeN;
        chargeNP = chargeN * chargeP;
        
        chargePP = chargeP * chargeP;        
	}

    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair){
		double sum = 0.0;
		double r2 = 0.0;

		IMolecule nitrogena = pair.getMolecule(0);
		IMolecule nitrogenb = pair.getMolecule(1);
		
		// to compute the midpoint distance between the two
		IVectorMutable com1 = (nitrogena.getChildList().getAtom(2)).getPosition();
		IVectorMutable com2 = (nitrogenb.getChildList().getAtom(2)).getPosition();
		
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
		
//		if (r2 > rC*rC){ 
//			return 0.0;
//		}
		//if(r2<1.6) return Double.POSITIVE_INFINITY;
		
		/*
		 * for the point/ atomic assignment
		 * refer to SpeciesN2ShellModel.java class
		 * 
		 */
        IVectorMutable Pan1 = nitrogena.getChildList().getAtom(0).getPosition();                                                                        
        IVectorMutable Pan2 = nitrogena.getChildList().getAtom(1).getPosition();
        IVectorMutable Pac  = nitrogena.getChildList().getAtom(2).getPosition();
        IVectorMutable Pap1 = nitrogena.getChildList().getAtom(3).getPosition();                                                                        
        IVectorMutable Pap2 = nitrogena.getChildList().getAtom(4).getPosition();

        IVectorMutable Pbn1 = nitrogenb.getChildList().getAtom(0).getPosition();
        IVectorMutable Pbn2 = nitrogenb.getChildList().getAtom(1).getPosition();
        IVectorMutable Pbc  = nitrogenb.getChildList().getAtom(2).getPosition();
        IVectorMutable Pbp1 = nitrogenb.getChildList().getAtom(3).getPosition();
        IVectorMutable Pbp2 = nitrogenb.getChildList().getAtom(4).getPosition();

        
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
    		 *   there is no interaction between the (-ve)charges 
    		 */
    		
    		for (int i=2; i<5; i++){
    			IVectorMutable dist = (nitrogena.getChildList().getAtom(i)).getPosition();
    			
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
    			IVectorMutable dist = (nitrogenb.getChildList().getAtom(i)).getPosition();
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
    
    private double calcDisperOverlap(double r2, double alpha, double epsilon, double delta){
    	double r = Math.sqrt(r2);
    	double a = 6.0/alpha;
    	double b = r/delta;
    	double b2 = b*b;
    	double b6 = b2*b2*b2;
    	
    	return (epsilon/(1-a))*(a*Math.exp(alpha*(1-b))-(1/b6));
    }
    
	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}
    

    private static final long serialVersionUID = 1L;
    
	protected IBoundary boundary;
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
	
	protected final IVectorMutable work, shift;

}
