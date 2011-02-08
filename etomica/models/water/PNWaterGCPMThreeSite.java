package etomica.models.water;

import etomica.api.IAtomList;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.math.SpecialFunctions;
import etomica.space.ISpace;

/**
 * GCPM Water potential with three association site for Wertheim approach.  
 * 
 * @author Hye Min Kim
 */
public class PNWaterGCPMThreeSite extends PNWaterGCPM {
	
	private double mincosThetaH;
	private int bondType;
	public static final int BONDAC = 1, BONDBC = 2, BONDCA = 3, BONDCB = 4;
	protected boolean isAssociation;
	double minE;

	public PNWaterGCPMThreeSite(ISpace space, double minE, boolean isAssociation) {
		super(space);
		this.isAssociation = isAssociation;
		this.minE = minE;
	}
	
	public double energy(IMoleculeList atoms){
        IAtomList water1Atoms = atoms.getMolecule(0).getChildList();
        IAtomList water2Atoms = atoms.getMolecule(1).getChildList();

        IVectorMutable O1r = water1Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        IVectorMutable O2r = water2Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        
        work.Ev1Mv2(O1r, O2r);
        
        double r2 = work.squared();
        
        if(r2<=core || r2 >= 12.25) {//outer shell = 3.5A
        	if(isAssociation) return 0.0;
        	else return super.energy(atoms);
        } 
        
        work.normalize();

        IVectorMutable H11r = water1Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        IVectorMutable H12r = water1Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();
        IVectorMutable H21r = water2Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        IVectorMutable H22r = water2Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();

        IVectorMutable M1r = water1Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        IVectorMutable M2r = water2Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        
        IVectorMutable rH11O2 = space.makeVector();
		rH11O2.E(H11r);
		rH11O2.ME(O2r);
		double distanceH11O2 = Math.sqrt(rH11O2.squared());
		
		IVectorMutable rH12O2 = space.makeVector();
		rH12O2.E(H12r);
		rH12O2.ME(O2r);
		double distanceH12O2 = Math.sqrt(rH12O2.squared());
		
		IVectorMutable rO1H21 = space.makeVector();
		rO1H21.E(O1r);
		rO1H21.ME(H21r);
		double distanceO1H21 = Math.sqrt(rO1H21.squared());
		
		IVectorMutable rO1H22 = space.makeVector();
		rO1H22.E(O1r);
		rO1H22.ME(H22r);
		double distanceO1H22 = Math.sqrt(rO1H22.squared());
		
        double[] distance = {distanceH11O2, distanceH12O2, distanceO1H21, distanceO1H22};
        double minDistance = distance[0];
        int minBond = 0;
        
        for (int i = 1; i < distance.length; i++){
        	if (distance[i] < minDistance){
        		minDistance = distance[i];
        		minBond = i;
        	}
        }
        
		IVectorMutable O1H11 = space.makeVector();
		IVectorMutable O1H12 = space.makeVector();
		IVectorMutable O2H11 = space.makeVector();
		IVectorMutable O2H12 = space.makeVector();
		IVectorMutable O2H21 = space.makeVector();
		IVectorMutable O1H1A = space.makeVector();
		IVectorMutable O2H2A = space.makeVector();
		IVectorMutable O1H21 = space.makeVector();
		IVectorMutable O1H22 = space.makeVector();
		IVectorMutable O2H22 = space.makeVector();
		IVectorMutable O1O2 = space.makeVector();
		IVectorMutable O2O1 = space.makeVector();
		IVectorMutable O1OPrime = space.makeVector();
		IVectorMutable O2OPrime = space.makeVector();
		IVectorMutable O2H11Prime = space.makeVector();
		IVectorMutable O2H12Prime = space.makeVector();
		IVectorMutable O1H21Prime = space.makeVector();
		IVectorMutable O1H22Prime = space.makeVector();
		IVectorMutable crossProduct1 = space.makeVector();
		IVectorMutable crossProduct2 = space.makeVector();
		IVectorMutable projectionO1 = space.makeVector();
		IVectorMutable projectionO2 = space.makeVector();
		IVectorMutable projectionH11 = space.makeVector();
		IVectorMutable projectionH12 = space.makeVector();
		IVectorMutable projectionH21 = space.makeVector();
		IVectorMutable projectionH22 = space.makeVector();
		O1H11.E(H11r);
		O1H11.ME(O1r);
		O1H11.normalize();
		O1H12.E(H12r);
		O1H12.ME(O1r);
		O1H12.normalize();
		
		crossProduct1.E(O1H11);
		crossProduct1.XE(O1H12);
		crossProduct1.normalize();
		
		O1H1A.E(O1H11);
		O1H1A.PE(O1H12);
		O1H1A.TE(-1);
		O1H1A.normalize();
		
		O1H21.E(H21r);
		O1H21.ME(O1r);
		O1H21.normalize();
		O1H22.E(H22r);
		O1H22.ME(O1r);
		O1H22.normalize();
		O2H11.E(H11r);
		O2H11.ME(O2r);
		O2H11.normalize();
		O2H12.E(H12r);
		O2H12.ME(O2r);
		O2H12.normalize();
		O1O2.E(O2r);
		O1O2.ME(O1r);
		O1O2.normalize();
		O2O1.Ea1Tv1(-1, O1O2);
		O2H21.E(H21r);
		O2H21.ME(O2r);
		O2H21.normalize();
		O2H22.E(H22r);
		O2H22.ME(O2r);
		O2H22.normalize();
		
		crossProduct2.E(O2H21);
		crossProduct2.XE(O2H22);
		crossProduct2.normalize();
		
		O2H2A.E(O2H21);
		O2H2A.PE(O2H22);
		O2H2A.TE(-1);
		O2H2A.normalize();
		
		double fractionO1O2 = crossProduct1.dot(O1O2);
		double fractionO2O1 = crossProduct2.dot(O2O1);
		double fractionO1H21 = crossProduct1.dot(O1H21);
		double fractionO1H22 = crossProduct1.dot(O1H22);
		double fractionO2H11 = crossProduct2.dot(O2H11);
		double fractionO2H12 = crossProduct2.dot(O2H12);
		projectionO1.Ea1Tv1(fractionO2O1,crossProduct2); 
		projectionO2.Ea1Tv1(fractionO1O2, crossProduct1);
		projectionH11.Ea1Tv1(fractionO2H11, crossProduct2);
		projectionH12.Ea1Tv1(fractionO2H12, crossProduct2);
		projectionH21.Ea1Tv1(fractionO1H21, crossProduct1);
		projectionH22.Ea1Tv1(fractionO1H22, crossProduct1);
		O1OPrime.Ev1Mv2(O1O2, projectionO2);
		O2OPrime.Ev1Mv2(O2O1, projectionO1);
		O1H21Prime.Ev1Mv2(O1H21,projectionH21);
		O1H22Prime.Ev1Mv2(O1H22,projectionH22);
		O2H11Prime.Ev1Mv2(O2H11,projectionH11);
		O2H12Prime.Ev1Mv2(O2H12,projectionH12);
		
		double cosTheta1=0;
		double cosTheta2=0;
		double cosTheta3=0;
		double cosTheta4=0;
		
		if (minBond == 0){
			cosTheta1 = O1OPrime.dot(O1H11);
        	cosTheta2 = O1OPrime.dot(O1O2);
        	cosTheta3 = O2H11Prime.dot(O2H2A);
        	cosTheta4 = O2H11Prime.dot(O2H11);
		}
		if (minBond == 1){
			cosTheta1 = O1OPrime.dot(O1H12);
        	cosTheta2 = O1OPrime.dot(O1O2);
        	cosTheta3 = O2H12Prime.dot(O2H2A);
        	cosTheta4 = O2H12Prime.dot(O2H12);
		}
		if (minBond == 2){
			cosTheta1 = O2OPrime.dot(O2H21);
        	cosTheta2 = O2OPrime.dot(O2O1);
        	cosTheta3 = O1H21Prime.dot(O1H1A);
        	cosTheta4 = O1H21Prime.dot(O1H21);
		}
		if (minBond == 3){
			cosTheta1 = O2OPrime.dot(O2H22);
        	cosTheta2 = O2OPrime.dot(O2O1);
        	cosTheta3 = O1H22Prime.dot(O1H1A);
        	cosTheta4 = O1H22Prime.dot(O1H22);
		}
        
        double energy = super.energy(atoms); 
     
        if (isAssociation&&bondType !=(minBond+1)){
        	return 0.0;
        }
//        if (isAssociation && (cosTheta1<0.6||cosTheta2<0.4||cosTheta3<-0.5||cosTheta4<0.00046)){
//        	return 0.0;
//        }
        

        if (isAssociation){
        	if (energy > -minE){
        		return 0.0;
        	}
        }
        else if(energy < -minE){
			return 0.0;
		}
 
		return energy;
	}
        
        public double getTheta() {return Math.acos(mincosThetaH);}
        
        /**
         * Accessor method for angle (in radians) describing width of cone.
         */
        public void setTheta(double t) {
            mincosThetaH   = Math.cos(t);
        }
        
        public void setBondType(int a){
        	bondType = a;
        }
}
