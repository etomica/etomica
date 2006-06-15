/*
 * Created on May 8, 2006
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.meam;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.potential.PotentialN;
import etomica.potential.PotentialSoft;
import etomica.space.NearestImageTransformer;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;

/**
 * @author ub2092
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class PotentialMEAM extends PotentialN implements PotentialSoft {
	
	public PotentialMEAM(Space space, ParameterSetMEAM p) {
		super(space);
        this.p = p;
        gradEi[0] = (Vector3D)space.makeVector();
    }

	/* (non-Javadoc)
	 * @see etomica.potential.Potential#getRange()
	 */
	public double getRange() {
		//from MEAMP2
		return Double.POSITIVE_INFINITY;
	}

	/* (non-Javadoc)
	 * @see etomica.potential.Potential#energy(etomica.atom.AtomSet)
	 */
	public double energy(AtomSet atoms) {
		AtomLeaf atom0 = (AtomLeaf)atoms.getAtom(0);
		
		//System.out.println(atom0);
        
        double sumRhoj0 = 0, sumRhoj1 = 0, sumRhoj2 = 0, sumRhoj3 = 0, 
			sumt1Rhoj0 = 0, sumt2Rhoj0 = 0, sumt3Rhoj0 = 0, 
			sumRhoj1x = 0, sumRhoj1y = 0, sumRhoj1z = 0,  
			sumRhoj2xx = 0, sumRhoj2xy = 0, sumRhoj2xz = 0,
			sumRhoj2yy = 0, sumRhoj2yz = 0, sumRhoj2zz = 0,
			sumRhoj3xxx = 0, sumRhoj3xxy = 0, sumRhoj3xxz = 0, 
			sumRhoj3xyy = 0, sumRhoj3xyz = 0, sumRhoj3xzz = 0,
			sumRhoj3yyy = 0, sumRhoj3yyz = 0, sumRhoj3yzz = 0, 
			sumRhoj3zzz = 0,
			sumPhi = 0;
		
        
        for(int j = 1; j < atoms.count(); j++) {
        	AtomLeaf atomj = (AtomLeaf) atoms.getAtom(j);
            rij.Ev1Mv2(atomj.coord.position(), atom0.coord.position());
            nearestImageTransformer.nearestImage(rij);
            double r = Math.sqrt(rij.squared());
            
            if (r > 4.0) continue; //Sn
            //if (r > 3.0) continue; //Cu
            
            //To determine amount of screening between atoms i and j 
            //by any atom k which may be between them
            double Sij = 1;
            for(int k = 1; k < atoms.count(); k++) {
            	AtomLeaf atomk = (AtomLeaf) atoms.getAtom(k);
            	if (k == j) continue;
            	rik.Ev1Mv2(atomk.coord.position(), atom0.coord.position());
            	nearestImageTransformer.nearestImage(rik);
            	double ik = Math.sqrt(rik.squared());
            	if (ik > 4.0*1.14) continue;
            	
            	double anglekij = Math.toDegrees(Math.acos(
            			((rij.x(0)*rik.x(0)) + 
            			 (rij.x(1)*rik.x(1)) + 
						 (rij.x(2)*rik.x(2)))
						 /(r*ik)));
            	
            	if (anglekij >= 90) continue;
            	
            	rkj.Ev1Mv2(atomk.coord.position(), atomj.coord.position());
            	nearestImageTransformer.nearestImage(rkj);
            	double kj = Math.sqrt(rkj.squared());
            	
            	//from Baskes, Angelo, & Bisson (1994)
            	double xik = (ik/kj)*(ik/kj);
            	double xjk = (kj/r)*(kj/r);
            	double C = ((2*(xik + xjk)) - ((xik - xjk)*(xik - xjk))- 1)/
						   (1 - ((xik - xjk)*(xik - xjk)));
            	
            	double Sijk;
            	if (C <= p.Cmin) { 
            		Sij = 0;
            		break;
            	}
            	else if (C >= p.Cmax) { //Sikj = 1, value of Sij won't change
            			continue;
            		}
            	else {
            			Sijk = Math.exp(-(((p.Cmax - C)/(C - p.Cmin))
            					*((p.Cmax - C)/(C - p.Cmin))));
            	}
            	
				Sij *= Sijk;
			
            }
    	
        if (Sij == 0) continue;    
       
    	//
    	//These rhoj terms are required for both the many-body and the 
    	//pair potential.  They are the contributions from each neighbor j
    	//to the partial electron densities of atom i.  rhoj0 is the contribution
    	//to the electron density in the s orbitals; rhoj1 is that in the p orbitals;
    	// rhoj2 is that in the d orbitals; rhoj3 is that in the f orbitals. 
    	//
    	//When the atoms in the pair being examined are different elements,
    	//the rhoj terms are not equivalent for both, because the constants
    	//for the elements in the expression for rhoj differ. This class functions
    	//for unmixed pairs of atoms only.
    	//
    	double rhoj0 = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0)) * Sij;
    	sumRhoj0 += (p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0))) * Sij;
    	double rhoj1 = p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0)) * Sij;
        sumRhoj1 += (p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0))) * Sij;
        double rhoj2 = p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0)) * Sij;
        sumRhoj2 += (p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0))) * Sij;
    	double rhoj3 = p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0)) * Sij;
    	sumRhoj3 += (p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0))) * Sij;
    	
    	sumt1Rhoj0 += (rhoj0 * p.t1);
    	sumt2Rhoj0 += (rhoj0 * p.t2);
    	sumt3Rhoj0 += (rhoj0 * p.t3);
    	
    	unitVector.E(rij);
        unitVector.normalize();
        
        double x = unitVector.x(0);
        double y = unitVector.x(1);
        double z = unitVector.x(2);
       
        sumRhoj1x += rhoj1 * x;
        sumRhoj1y += rhoj1 * y;
        sumRhoj1z += rhoj1 * z;  	
        
        sumRhoj2xx += rhoj2 * x * x;
        sumRhoj2xy += rhoj2 * x * y;
        sumRhoj2xz += rhoj2 * x * z;
        sumRhoj2yy += rhoj2 * y * y;
        sumRhoj2yz += rhoj2 * y * z;
        sumRhoj2zz += rhoj2 * z * z;
        
        sumRhoj3xxx += rhoj3 * x * x * x;
        sumRhoj3xxy += rhoj3 * x * x * y;
        sumRhoj3xxz += rhoj3 * x * x * z;
        sumRhoj3xyy += rhoj3 * x * y * y;
        sumRhoj3xyz += rhoj3 * x * y * z;
        sumRhoj3xzz += rhoj3 * x * z * z;
        sumRhoj3yyy += rhoj3 * y * y * y;
        sumRhoj3yyz += rhoj3 * y * y * z;
        sumRhoj3yzz += rhoj3 * y * z * z;
        sumRhoj3zzz += rhoj3 * z * z * z;
        
        //to calculate the repulsive pair potential, phi
        //Should rhoj0 within phi contain Sij? Phi itself is multiplied by Sij...
        double rhoj0Ref = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0));
        double rhoiRef = p.Z * rhoj0Ref;   //FCC reference structure, Z = 12
        double a = p.alpha * ((r/p.r0) - 1.0);
    	double EuRef = - p.Ec * (1.0 + a) * Math.exp(-a);
    	double FRef = p.A * p.Ec * (rhoiRef/p.Z) * Math.log(rhoiRef/p.Z);
    	sumPhi += ((2.0/p.Z) * (EuRef - FRef)) * Sij;
        }
        
        double rhoi0 = sumRhoj0;
        
        double rhoi1 = Math.sqrt( (sumRhoj1x * sumRhoj1x)
    				+ (sumRhoj1y * sumRhoj1y)
					+ (sumRhoj1z * sumRhoj1z));
       
        double rhoi2 = Math.sqrt((sumRhoj2xx * sumRhoj2xx)
    			+ (2.0 * sumRhoj2xy * sumRhoj2xy)
				+ (2.0 * sumRhoj2xz * sumRhoj2xz)
    			+ (sumRhoj2yy * sumRhoj2yy) 
				+ (2.0 * sumRhoj2yz * sumRhoj2yz)
    			+ (sumRhoj2zz * sumRhoj2zz) 
				- ((1.0/3.0) * sumRhoj2 * sumRhoj2));
        
        double rhoi3 = Math.sqrt( 
        		   (sumRhoj3xxx * sumRhoj3xxx)
    		+(3.0 * sumRhoj3xxy * sumRhoj3xxy)
    		+(3.0 * sumRhoj3xxz * sumRhoj3xxz)
    		+(3.0 * sumRhoj3xyy * sumRhoj3xyy)
    		+(6.0 * sumRhoj3xyz * sumRhoj3xyz)
    		+(3.0 * sumRhoj3xzz * sumRhoj3xzz)
    		+      (sumRhoj3yyy * sumRhoj3yyy)
    		+(3.0 * sumRhoj3yyz * sumRhoj3yyz)
    		+(3.0 * sumRhoj3yzz * sumRhoj3yzz)
    		+      (sumRhoj3zzz * sumRhoj3zzz));
        
        double tav1 = sumt1Rhoj0 / sumRhoj0;

        double tav2 = sumt2Rhoj0 / sumRhoj0;
        
        double tav3 = sumt3Rhoj0 / sumRhoj0;
        
        double gamma = (tav1 * (rhoi1/rhoi0) * (rhoi1/rhoi0))
    		+ (tav2 * (rhoi2/rhoi0) * (rhoi2/rhoi0))
    		+ (tav3 * (rhoi3/rhoi0) * (rhoi3/rhoi0));

        //The following expression for the background electron density of atom i
		//is appropriate for Sn only.
        double rhoi = (2.0 * rhoi0) / (1.0 + Math.exp(-gamma)); //Sn
        //double rhoi = rhoi0 * Math.sqrt(1.0 + gamma); //Cu
        
		double F = p.A * p.Ec * (rhoi/p.Z) * Math.log(rhoi/p.Z);
		
		return F + (0.5*sumPhi);
	}

	/* (non-Javadoc)
	 * @see etomica.potential.Potential#setPhase(etomica.phase.Phase)
	 */
	public void setPhase(Phase phase) {
		nearestImageTransformer = phase.getBoundary();

	}

	/* (non-Javadoc)
	 * @see etomica.potential.PotentialSoft#virial(etomica.atom.AtomSet)
	 */
	public double virial(AtomSet atoms) {
		return 0.0;
	}

	/* (non-Javadoc)
	 * @see etomica.potential.PotentialSoft#gradient(etomica.atom.AtomSet)
	 */
	public Vector[] gradient(AtomSet atoms) {
        AtomLeaf atom0 = (AtomLeaf)atoms.getAtom(0);
        
        double sumRhoj0 = 0, sumRhoj1 = 0, sumRhoj2 = 0, sumRhoj3 = 0, 
		sumt1Rhoj0 = 0, sumt2Rhoj0 = 0, sumt3Rhoj0 = 0, 
		sumRhoj1x = 0, sumRhoj1y = 0, sumRhoj1z = 0,  
		sumRhoj2xx = 0, sumRhoj2xy = 0,  sumRhoj2xz = 0,
		sumRhoj2yy = 0, sumRhoj2yz = 0, sumRhoj2zz = 0,
		sumRhoj3xxx = 0, sumRhoj3xxy = 0, sumRhoj3xxz = 0, 
		sumRhoj3xyy = 0, sumRhoj3xyz = 0, sumRhoj3xzz = 0,
		sumRhoj3yyy = 0, sumRhoj3yyz = 0, sumRhoj3yzz = 0, sumRhoj3zzz = 0,
		sumPhi = 0; //sumPhi may not be needed in gradient() method
        
        sumGradRhoj0.E(0);
        sumGradRhoj1.E(0);
        sumGradRhoj2.E(0);
        sumGradRhoj3.E(0);
        sumGradRhoj1x.E(0);
        sumGradRhoj1y.E(0);
        sumGradRhoj1z.E(0);
        sumGradRhoj2xx.E(0);
        sumGradRhoj2xy.E(0);
        sumGradRhoj2xz.E(0);
        sumGradRhoj2yy.E(0);
        sumGradRhoj2yz.E(0);
        sumGradRhoj2zz.E(0);
        sumGradRhoj3xxx.E(0);
        sumGradRhoj3xxy.E(0);
        sumGradRhoj3xxz.E(0);
        sumGradRhoj3xyy.E(0);
        sumGradRhoj3xyz.E(0);
        sumGradRhoj3xzz.E(0);
        sumGradRhoj3yyy.E(0);
        sumGradRhoj3yyz.E(0);
        sumGradRhoj3yzz.E(0);
        sumGradRhoj3zzz.E(0);
        sumt1GradRhoj0.E(0);
        sumt2GradRhoj0.E(0);
        sumt3GradRhoj0.E(0);
        sumGradPhi.E(0);
        
        for(int j = 1; j < atoms.count(); j++) {
        	
        	AtomLeaf atomj = (AtomLeaf) atoms.getAtom(j);
            rij.Ev1Mv2(atomj.coord.position(), atom0.coord.position());
            nearestImageTransformer.nearestImage(rij);
            double r = Math.sqrt(rij.squared());
            
            if (r > 4.0) continue; //Sn
            //if (r > 3.0) continue; //Cu
            
            //To determine amount of screening between atoms i and j 
            //by any atom k which may be between them
            double Sij = 1;
            giSij.E(0);
            for(int k = 1; k < atoms.count(); k++) {
            	AtomLeaf atomk = (AtomLeaf) atoms.getAtom(k);
            	if (k == j) continue;
            	rik.Ev1Mv2(atomk.coord.position(), atom0.coord.position());
            	nearestImageTransformer.nearestImage(rik);
            	double ik = Math.sqrt(rik.squared());
            	if (ik > 4.0*1.14) continue;
            	
            	double anglekij = Math.toDegrees(Math.acos(
            			((rij.x(0)*rik.x(0)) + 
            			 (rij.x(1)*rik.x(1)) + 
						 (rij.x(2)*rik.x(2)))
						 /(r*ik)));
            	
            	if (anglekij >= 90) continue;
            	
            	rkj.Ev1Mv2(atomk.coord.position(), atomj.coord.position());
            	nearestImageTransformer.nearestImage(rkj);
            	double kj = Math.sqrt(rkj.squared());
            	
            	//from Baskes, Angelo, & Bisson (1994)
            	double xik = (ik/kj)*(ik/kj);
            	double xkj = (kj/r)*(kj/r);
            	
            	double C = ( (2*(xik + xkj)) - ((xik - xkj)*(xik - xkj))- 1 )
							/ (1 - ((xik - xkj)*(xik - xkj)));
            	
            	double Sijk;
            	if (C <= p.Cmin) { 
            		Sij = 0;
            		break;
            	}
            	else if (C >= p.Cmax) { //Sikj = 1, value of Sij won't change
            			continue;
            		}
            	else {
            			Sijk = Math.exp(-(((p.Cmax - C)/(C - p.Cmin))
            					*((p.Cmax - C)/(C - p.Cmin))));
            	}
            	
            	giRij.Ea1Tv1(-1/r, rij);
            	gjRij.Ea1Tv1(-1, giRij);
            	gkRij.E(0);
            	
            	giRik.Ea1Tv1(-1/ik, rik);
            	gjRik.E(0);
            	gkRik.Ea1Tv1(-1, giRik);
            	
            	giRkj.E(0);
            	gjRkj.Ea1Tv1(-1/kj, rkj);
            	gkRkj.Ea1Tv1(-1, gjRkj);
            	
            	giXik.Ea1Tv1(-ik/(kj*kj), giRkj);
            	giXik.PEa1Tv1(1/kj, giRik);
            	giXik.TE(2*ik/kj);
            	
            	giXkj.Ea1Tv1(-kj/(r*r), giRij);
            	giXkj.PEa1Tv1(1/r, giRkj);
            	giXkj.TE(2*kj/r);
            	
            	gjXik.Ea1Tv1(-ik/(kj*kj), gjRkj);
            	gjXik.PEa1Tv1(1/kj, gjRik);
            	gjXik.TE(2*ik/kj);
            	
            	gjXkj.Ea1Tv1(-kj/(r*r), gjRij);
            	gjXkj.PEa1Tv1(1/r, gjRkj);
            	gjXkj.TE(2*kj/r);
            	
            	gkXik.Ea1Tv1(-ik/(kj*kj), gkRkj);
            	gkXik.PEa1Tv1(1/kj, gkRik);
            	gkXik.TE(2*ik/kj);
            	
            	gkXkj.Ea1Tv1(-kj/(r*r), gkRij);
            	gkXkj.PEa1Tv1(1/r, gkRkj);
            	gkXkj.TE(2*kj/r);
            	
				giC.Ea1Tv1( 1 + (xik - xkj)*(C - 1), giXik);
		    	giC.PEa1Tv1(1 - (xik - xkj)*(C + 1), giXkj);
		    	giC.TE( 2 / ( 1 - ((xik - xkj)*(xik - xkj)) ));
		    	
		    	gjC.Ea1Tv1( 1 + (xik - xkj)*(C - 1), gjXik);
		    	gjC.PEa1Tv1(1 - (xik - xkj)*(C + 1), gjXkj);
		    	gjC.TE( 2 / ( 1 - ((xik - xkj)*(xik - xkj)) ));
		    	
		    	gkC.Ea1Tv1( 1 + (xik - xkj)*(C - 1), gkXik);
		    	gkC.PEa1Tv1(1 - (xik - xkj)*(C + 1), gkXkj);
		    	gkC.TE( 2 / ( 1 - ((xik - xkj)*(xik - xkj)) ));
		    	
		    	giSijk.Ea1Tv1( 2*Sijk*(p.Cmax - C)/((C - p.Cmin)*(C-p.Cmin) )
		    			* ( ((p.Cmax - C)/(C - p.Cmin)) + 1 ), giC);
		    	
		    	//The Sij value used to calculate gradSij is that for previous k's, 
		    	//or, for the first k considered, 1.  Same goes for the gradSij value
		    	//itself, except it's initialized as the zero vector...
		    	giSij.TE(Sijk);
		    	giSij.PEa1Tv1(Sij, giSijk);
		    
		    	Sij *= Sijk;
		    	
            }
    	
            if (Sij == 0) continue;
        	
        	double rhoj0 = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0)) * Sij;
        	sumRhoj0 += p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0)) * Sij;
        	double rhoj1 = p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0)) * Sij;
            sumRhoj1 += p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0)) * Sij;
            double rhoj2 = p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0)) * Sij;
            sumRhoj2 += p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0)) * Sij;
        	double rhoj3 = p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0)) * Sij;
        	sumRhoj3 += p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0)) * Sij;
        	
        	sumt1Rhoj0 += rhoj0 * p.t1;
        	sumt2Rhoj0 += rhoj0 * p.t2;
        	sumt3Rhoj0 += rhoj0 * p.t3;
        	
        	unitVector.E(rij);
            unitVector.normalize();
            
            double x = unitVector.x(0);
            double y = unitVector.x(1);
            double z = unitVector.x(2);
           
            sumRhoj1x += rhoj1 * x;
            sumRhoj1y += rhoj1 * y;
            sumRhoj1z += rhoj1 * z;  	
            
            sumRhoj2xx += rhoj2 * x * x;
            sumRhoj2xy += rhoj2 * x * y;
            sumRhoj2xz += rhoj2 * x * z;
            sumRhoj2yy += rhoj2 * y * y;
            sumRhoj2yz += rhoj2 * y * z;
            sumRhoj2zz += rhoj2 * z * z;
            
            sumRhoj3xxx += rhoj3 * x * x * x;
            sumRhoj3xxy += rhoj3 * x * x * y;
            sumRhoj3xxz += rhoj3 * x * x * z;
            sumRhoj3xyy += rhoj3 * x * y * y;
            sumRhoj3xyz += rhoj3 * x * y * z;
            sumRhoj3xzz += rhoj3 * x * z * z;
            sumRhoj3yyy += rhoj3 * y * y * y;
            sumRhoj3yyz += rhoj3 * y * y * z;
            sumRhoj3yzz += rhoj3 * y * z * z;
            sumRhoj3zzz += rhoj3 * z * z * z;
        
            //to calculate the repulsive pair potential, phi
            //Should rhoj0 within phi contain Sij? Phi itself is multiplied by Sij...
            //FCC reference structure, Z = 12
            double rhoj0Ref = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0));
            double rhoiRef = p.Z * rhoj0Ref;   
            double a = p.alpha * ((r/p.r0) - 1.0);
        	double EuRef = - p.Ec * (1.0 + a) * Math.exp(-a);
        	double FRef = p.A * p.Ec * (rhoiRef/p.Z) * Math.log(rhoiRef/p.Z);
        	double phi = ((2.0/p.Z) * (EuRef - FRef)) * Sij;
        	sumPhi += ((2.0/p.Z) * (EuRef - FRef)) * Sij;
        	
    	
    	vector100.setX(0,1.0);
    	vector010.setX(1,1.0);
    	vector001.setX(2,1.0);
    	
    	gradx.Ea1Tv1(x/(r*r),rij);
    	gradx.PEa1Tv1(-1.0/r, vector100);
    	
    	grady.Ea1Tv1(y/(r*r),rij);
    	grady.PEa1Tv1(-1.0/r, vector010);
    	
    	gradz.Ea1Tv1(z/(r*r),rij);
    	gradz.PEa1Tv1(-1.0/r, vector001);
    	
    	gradr.Ea1Tv1(-1.0/r,rij);
    	
    	//Calculations to determine gradient of phi
    	//FCC reference structure, Z = 12
    	gradRhoj0Ref.Ea1Tv1(-rhoj0Ref*p.beta0/p.r0, gradr);
    	gradRhoiRef.Ea1Tv1(p.Z, gradRhoj0Ref);
    	gradFRef.Ea1Tv1( (p.A*p.Ec/p.Z)*
			(1.0 + Math.log(rhoiRef/p.Z)), gradRhoiRef);
    	gradERef.Ea1Tv1( (p.Ec*p.alpha*p.alpha/p.r0) * ((r/p.r0) - 1.0)
					*(Math.exp(-p.alpha*((r/p.r0)-1.0))), gradr);
	
    	gradPhi.E(gradERef);
    	gradPhi.ME(gradFRef);
    	gradPhi.TE(2.0 * Sij /p.Z);
    	gradPhi.PEa1Tv1(phi/Sij, giSij);
    	sumGradPhi.PE(gradPhi);
    
    	//Gradient of rhoj0
    	gradRhoj0.Ea1Tv1(rhoj0/Sij, giSij);
    	gradRhoj0.PEa1Tv1(-rhoj0*p.beta0/(p.r0), gradr);
    	sumGradRhoj0.PE(gradRhoj0);
    	
    	//Gradient of rhoj1
    	gradRhoj1.Ea1Tv1(rhoj1/Sij, giSij);
    	gradRhoj1.PEa1Tv1(-rhoj1*p.beta1/(p.r0), gradr);
    	sumGradRhoj1.PE(gradRhoj1);
    	
    	//Gradient of rhoj2
    	gradRhoj2.Ea1Tv1(rhoj2/Sij, giSij);
    	gradRhoj2.PEa1Tv1(-rhoj2*p.beta2/(p.r0), gradr);
    	sumGradRhoj2.PE(gradRhoj2);
    	
    	//Gradient of rhoj3
    	gradRhoj3.Ea1Tv1(rhoj3/Sij, giSij);
    	gradRhoj3.PEa1Tv1(-rhoj3*p.beta3/(p.r0), gradr);
    	sumGradRhoj3.PE(gradRhoj3);
    	
    	//Gradient of rhoj1x
    	gradRhoj1x.Ea1Tv1(rhoj1, gradx);
    	gradRhoj1x.PEa1Tv1(x, gradRhoj1);
    	sumGradRhoj1x.PE(gradRhoj1x);
    	
    	//Gradient of rhoj1y
    	gradRhoj1y.Ea1Tv1(rhoj1, grady);
    	gradRhoj1y.PEa1Tv1(y, gradRhoj1);
    	sumGradRhoj1y.PE(gradRhoj1y);
    	
    	//Gradient of rhoj1z
    	gradRhoj1z.Ea1Tv1(rhoj1, gradz);
    	gradRhoj1z.PEa1Tv1(z, gradRhoj1);
    	sumGradRhoj1z.PE(gradRhoj1z);
    	
    	//Gradient of rhoj2xx
    	gradRhoj2xx.Ea1Tv1(2.0*rhoj2*x, gradx);
    	gradRhoj2xx.PEa1Tv1(x*x, gradRhoj2);
    	sumGradRhoj2xx.PE(gradRhoj2xx);
    	
    	//Gradient of rhoj2xy
    	gradRhoj2xy.Ea1Tv1(rhoj2*x, grady);
    	gradRhoj2xy.PEa1Tv1(rhoj2*y, gradx);
    	gradRhoj2xy.PEa1Tv1(x*y, gradRhoj2);
    	sumGradRhoj2xy.PE(gradRhoj2xy);
    	
    	//Gradient of rhoj2xz
    	gradRhoj2xz.Ea1Tv1(rhoj2*x, gradz);
    	gradRhoj2xz.PEa1Tv1(rhoj2*z, gradx);
    	gradRhoj2xz.PEa1Tv1(x*z, gradRhoj2);
    	sumGradRhoj2xz.PE(gradRhoj2xz);
    	
    	//Gradient of rhoj2yy
    	gradRhoj2yy.Ea1Tv1(2.0*rhoj2*y, grady);
    	gradRhoj2yy.PEa1Tv1(y*y, gradRhoj2);
    	sumGradRhoj2yy.PE(gradRhoj2yy);
    	
    	//Gradient of rhoj2yz
    	gradRhoj2yz.Ea1Tv1(rhoj2*y, gradz);
    	gradRhoj2yz.PEa1Tv1(rhoj2*z, grady);
    	gradRhoj2yz.PEa1Tv1(y*z, gradRhoj2);
    	sumGradRhoj2yz.PE(gradRhoj2yz);
    	
    	//Gradient of rhoj2zz
    	gradRhoj2zz.Ea1Tv1(2.0*rhoj2*z, gradz);
    	gradRhoj2zz.PEa1Tv1(z*z, gradRhoj2);
    	sumGradRhoj2zz.PE(gradRhoj2zz);
    	
    	//Gradient of rhoj3xxx
    	gradRhoj3xxx.Ea1Tv1(3.0*rhoj3*x*x, gradx);
    	gradRhoj3xxx.PEa1Tv1(x*x*x, gradRhoj3);
    	sumGradRhoj3xxx.PE(gradRhoj3xxx);
    	
    	//Gradient of rhoj3xxy
    	gradRhoj3xxy.Ea1Tv1(rhoj3*x*x, grady);
    	gradRhoj3xxy.PEa1Tv1(2.0*rhoj3*x*y, gradx);
    	gradRhoj3xxy.PEa1Tv1(x*x*y, gradRhoj3);
    	sumGradRhoj3xxy.PE(gradRhoj3xxy);
    	
    	//Gradient of rhoj3xxz
    	gradRhoj3xxz.Ea1Tv1(rhoj3*x*x, gradz);
    	gradRhoj3xxz.PEa1Tv1(2.0*rhoj3*x*z, gradx);
    	gradRhoj3xxz.PEa1Tv1(x*x*z, gradRhoj3);
    	sumGradRhoj3xxz.PE(gradRhoj3xxz);
    	
    	//Gradient of rhoj3xyy
    	gradRhoj3xyy.Ea1Tv1(2.0*rhoj3*x*y, grady);
    	gradRhoj3xyy.PEa1Tv1(rhoj3*y*y, gradx);
    	gradRhoj3xyy.PEa1Tv1(x*y*y, gradRhoj3);
    	sumGradRhoj3xyy.PE(gradRhoj3xyy);
    	
    	//Gradient of rhoj3xyz
    	gradRhoj3xyz.Ea1Tv1(rhoj3*x*y, gradz);
    	gradRhoj3xyz.PEa1Tv1(rhoj3*x*z, grady);
    	gradRhoj3xyz.PEa1Tv1(rhoj3*y*z, gradx);
    	gradRhoj3xyz.PEa1Tv1(x*y*z, gradRhoj3);
    	sumGradRhoj3xyz.PE(gradRhoj3xyz);
    	
    	//Gradient of rhoj3xzz
    	gradRhoj3xzz.Ea1Tv1(2.0*rhoj3*x*z, gradz);
    	gradRhoj3xzz.PEa1Tv1(rhoj3*z*z, gradx);
    	gradRhoj3xzz.PEa1Tv1(x*z*z, gradRhoj3);
    	sumGradRhoj3xzz.PE(gradRhoj3xzz);
    	
    	//Gradient of rhoj3yyy
    	gradRhoj3yyy.Ea1Tv1(3.0*rhoj3*y*y, grady);
    	gradRhoj3yyy.PEa1Tv1(y*y*y, gradRhoj3);
    	sumGradRhoj3yyy.PE(gradRhoj3yyy);
    	
    	//Gradient of rhoj3yyz
    	gradRhoj3yyz.Ea1Tv1(rhoj3*y*y, gradz);
    	gradRhoj3yyz.PEa1Tv1(2.0*rhoj3*y*z, grady);
    	gradRhoj3yyz.PEa1Tv1(y*y*z, gradRhoj3);
    	sumGradRhoj3yyz.PE(gradRhoj3yyz);
    	
    	//Gradient of rhoj3yzz
    	gradRhoj3yzz.Ea1Tv1(2.0*rhoj3*y*z, gradz);
    	gradRhoj3yzz.PEa1Tv1(rhoj3*z*z, grady);
    	gradRhoj3yzz.PEa1Tv1(y*z*z, gradRhoj3);
    	sumGradRhoj3yzz.PE(gradRhoj3yzz);
    	
    	//Gradient of rhoj3zzz
    	gradRhoj3zzz.Ea1Tv1(3.0*rhoj3*z*z, gradz);
    	gradRhoj3zzz.PEa1Tv1(z*z*z, gradRhoj3);
    	sumGradRhoj3zzz.PE(gradRhoj3zzz);
    	
    	//t1GradRhoj0
    	t1GradRhoj0.Ea1Tv1(p.t1, gradRhoj0);
    	sumt1GradRhoj0.PE(t1GradRhoj0);
    	
    	//t2GradRhoj0
    	t2GradRhoj0.Ea1Tv1(p.t2, gradRhoj0);
    	sumt2GradRhoj0.PE(t2GradRhoj0);
    	
    	//t3GradRhoj0
    	t3GradRhoj0.Ea1Tv1(p.t3, gradRhoj0);
    	sumt3GradRhoj0.PE(t3GradRhoj0);
        }
        
        double rhoi0 = sumRhoj0;
        
        double rhoi1 = Math.sqrt( (sumRhoj1x * sumRhoj1x)
    				+ (sumRhoj1y * sumRhoj1y)+ (sumRhoj1z * sumRhoj1z));
       
        double rhoi2 = Math.sqrt((sumRhoj2xx * sumRhoj2xx)
    			+ (2.0 * sumRhoj2xy * sumRhoj2xy)+ (2.0 * sumRhoj2xz * sumRhoj2xz)
    			+ (sumRhoj2yy * sumRhoj2yy) + (2.0 * sumRhoj2yz * sumRhoj2yz)
    			+ (sumRhoj2zz * sumRhoj2zz) - ((1.0/3.0) * sumRhoj2 * sumRhoj2));
        
        double rhoi3 = Math.sqrt( 
        		   (sumRhoj3xxx * sumRhoj3xxx)
    		+(3.0 * sumRhoj3xxy * sumRhoj3xxy)
    		+(3.0 * sumRhoj3xxz * sumRhoj3xxz)
    		+(3.0 * sumRhoj3xyy * sumRhoj3xyy)
    		+(6.0 * sumRhoj3xyz * sumRhoj3xyz)
    		+(3.0 * sumRhoj3xzz * sumRhoj3xzz)
    		+      (sumRhoj3yyy * sumRhoj3yyy)
    		+(3.0 * sumRhoj3yyz * sumRhoj3yyz)
    		+(3.0 * sumRhoj3yzz * sumRhoj3yzz)
    		+      (sumRhoj3zzz * sumRhoj3zzz));
        
        double tav1 = sumt1Rhoj0 / sumRhoj0;

        double tav2 = sumt2Rhoj0 / sumRhoj0;
        
        double tav3 = sumt3Rhoj0 / sumRhoj0;
        
        double gamma = (tav1 * (rhoi1/rhoi0) * (rhoi1/rhoi0))
    		+ (tav2 * (rhoi2/rhoi0) * (rhoi2/rhoi0))
    		+ (tav3 * (rhoi3/rhoi0) * (rhoi3/rhoi0));

        //The following expression for the background electron density of atom i
		//is appropriate for Sn only.
        double rhoi = (2.0 * rhoi0) / (1.0 + Math.exp(-gamma)); //Sn
        
        //double rhoi = rhoi0 * Math.sqrt(1.0 + gamma); //Cu
        
        gradRhoi0.E(sumGradRhoj0);
		
		gradRhoi1.Ea1Tv1(sumRhoj1x, sumGradRhoj1x);
		gradRhoi1.PEa1Tv1(sumRhoj1y, sumGradRhoj1y);
		gradRhoi1.PEa1Tv1(sumRhoj1z, sumGradRhoj1z);
		gradRhoi1.TE(1.0/rhoi1);
		
		gradRhoi2.Ea1Tv1(sumRhoj2xx, sumGradRhoj2xx);
		gradRhoi2.PEa1Tv1(2.0 * sumRhoj2xy, sumGradRhoj2xy);
		gradRhoi2.PEa1Tv1(2.0 * sumRhoj2xz, sumGradRhoj2xz);
		gradRhoi2.PEa1Tv1(sumRhoj2yy, sumGradRhoj2yy);
		gradRhoi2.PEa1Tv1(2.0 * sumRhoj2yz, sumGradRhoj2yz);
		gradRhoi2.PEa1Tv1(sumRhoj2zz, sumGradRhoj2zz);
		gradRhoi2.PEa1Tv1(-(1.0/3.0) * sumRhoj2, sumGradRhoj2);
		gradRhoi2.TE(1.0/rhoi2);
		
		gradRhoi3.Ea1Tv1(sumRhoj3xxx, sumGradRhoj3xxx);
		gradRhoi3.PEa1Tv1(3.0 * sumRhoj3xxy, sumGradRhoj3xxy);
		gradRhoi3.PEa1Tv1(3.0 * sumRhoj3xxz, sumGradRhoj3xxz);
		gradRhoi3.PEa1Tv1(3.0 * sumRhoj3xyy, sumGradRhoj3xyy);
		gradRhoi3.PEa1Tv1(6.0 * sumRhoj3xyz, sumGradRhoj3xyz);
		gradRhoi3.PEa1Tv1(3.0 * sumRhoj3xzz, sumGradRhoj3xzz);
		gradRhoi3.PEa1Tv1(sumRhoj3yyy, sumGradRhoj3yyy);
		gradRhoi3.PEa1Tv1(3.0 * sumRhoj3yyz, sumGradRhoj3yyz);
		gradRhoi3.PEa1Tv1(3.0 * sumRhoj3yzz, sumGradRhoj3yzz);
		gradRhoi3.PEa1Tv1(sumRhoj3zzz, sumGradRhoj3zzz);
		gradRhoi3.TE(1.0/rhoi3);
		
		gradtav1.Ea1Tv1(-sumt1Rhoj0/(rhoi0*rhoi0), sumGradRhoj0);
		gradtav1.PEa1Tv1(1.0/rhoi0, sumt1GradRhoj0);
		
		gradtav2.Ea1Tv1(-sumt2Rhoj0/(rhoi0*rhoi0), sumGradRhoj0);
		gradtav2.PEa1Tv1(1.0/rhoi0, sumt2GradRhoj0);
		
		gradtav3.Ea1Tv1(-sumt3Rhoj0/(rhoi0*rhoi0), sumGradRhoj0);
		gradtav3.PEa1Tv1(1.0/rhoi0, sumt3GradRhoj0);
		
		gradGamma.Ea1Tv1( (-1.0/rhoi0)*( (tav1*rhoi1*rhoi1) + (tav2*rhoi2*rhoi2)
				+ (tav3*rhoi3*rhoi3) ), gradRhoi0);
		gradGamma.PEa1Tv1(tav1*rhoi1, gradRhoi1);
		gradGamma.PEa1Tv1(tav2*rhoi2, gradRhoi2);
		gradGamma.PEa1Tv1(tav3*rhoi3, gradRhoi3);
		gradGamma.TE(2.0);
		gradGamma.PEa1Tv1(rhoi1*rhoi1, gradtav1);
		gradGamma.PEa1Tv1(rhoi2*rhoi2, gradtav2);
		gradGamma.PEa1Tv1(rhoi3*rhoi3, gradtav3);
		gradGamma.TE(1.0/(rhoi0*rhoi0));
		
		//Sn
		
		gradRhoi.Ea1Tv1(rhoi0*Math.exp(-gamma)/(1.0 + Math.exp(-gamma)), gradGamma);
		gradRhoi.PE(gradRhoi0);
		gradRhoi.TE(2.0/(1.0+Math.exp(-gamma)));
		
		
		//Cu
		/**
		gradRhoi.Ea1Tv1(rhoi0*0.5*Math.sqrt(1.0/(1.0 + gamma)), gradGamma);
		gradRhoi.PEa1Tv1(Math.sqrt(1.0+gamma), gradRhoi0);
	    **/
		
		gradF.Ea1Tv1( (p.A*p.Ec/p.Z)*
				(1.0 + Math.log(rhoi/p.Z)), gradRhoi);
		
		//System.out.println("gradF is " + gradF);
		//System.exit(0);
		
		gradEi[0].E(gradF);
		gradEi[0].PEa1Tv1(0.5, sumGradPhi);
		
		return gradEi;
	}
	
	private final Vector3D unitVector = (Vector3D)space.makeVector();
    private final Vector3D oppositeDirection = (Vector3D)space.makeVector();
    private final Vector3D vector100 = (Vector3D)space.makeVector();
    private final Vector3D vector010 = (Vector3D)space.makeVector();
    private final Vector3D vector001 = (Vector3D)space.makeVector();
    private final Vector3D giRij = (Vector3D)space.makeVector();
    private final Vector3D gjRij = (Vector3D)space.makeVector();
    private final Vector3D gkRij = (Vector3D)space.makeVector();
    private final Vector3D giRik = (Vector3D)space.makeVector();
    private final Vector3D gjRik = (Vector3D)space.makeVector();
    private final Vector3D gkRik = (Vector3D)space.makeVector();
    private final Vector3D giRkj = (Vector3D)space.makeVector();
    private final Vector3D gjRkj = (Vector3D)space.makeVector();
    private final Vector3D gkRkj = (Vector3D)space.makeVector();
    private final Vector3D giXik = (Vector3D)space.makeVector();
    private final Vector3D gjXik = (Vector3D)space.makeVector();
    private final Vector3D gkXik = (Vector3D)space.makeVector();
    private final Vector3D giXkj = (Vector3D)space.makeVector();
    private final Vector3D gjXkj = (Vector3D)space.makeVector();
    private final Vector3D gkXkj = (Vector3D)space.makeVector();
    private final Vector3D giC = (Vector3D)space.makeVector();
    private final Vector3D gjC = (Vector3D)space.makeVector();
    private final Vector3D gkC = (Vector3D)space.makeVector();
    private final Vector3D giSijk = (Vector3D)space.makeVector();
    private final Vector3D giSij = (Vector3D)space.makeVector();
    
    private final Vector3D gradx = (Vector3D)space.makeVector();
    private final Vector3D grady = (Vector3D)space.makeVector();
    private final Vector3D gradz = (Vector3D)space.makeVector();
    private final Vector3D gradr = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj1 = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj1 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2 = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3 = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj1x = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj1x = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj1y = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj1y = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj1z = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj1z = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2xx = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2xx = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2xy = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2xy = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2xz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2xz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2yy = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2yy = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2yz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2yz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj2zz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj2zz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3xxx = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3xxx = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3xxy = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3xxy = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3xxz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3xxz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3xyy = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3xyy = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3xyz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3xyz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3xzz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3xzz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3yyy = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3yyy = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3yyz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3yyz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3yzz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3yzz = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj3zzz = (Vector3D)space.makeVector();
    private final Vector3D sumGradRhoj3zzz = (Vector3D)space.makeVector();
    private final Vector3D t1GradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D sumt1GradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D t2GradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D sumt2GradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D t3GradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D sumt3GradRhoj0 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoj0Ref = (Vector3D)space.makeVector();
    private final Vector3D gradRhoiRef = (Vector3D)space.makeVector();
    private final Vector3D gradERef = (Vector3D)space.makeVector();
    private final Vector3D gradFRef = (Vector3D)space.makeVector();
    private final Vector3D gradPhi = (Vector3D)space.makeVector();
    private final Vector3D sumGradPhi = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi0 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi1 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi2 = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi3 = (Vector3D)space.makeVector();
    private final Vector3D gradtav1 = (Vector3D)space.makeVector();
    private final Vector3D gradtav2 = (Vector3D)space.makeVector();
    private final Vector3D gradtav3 = (Vector3D)space.makeVector();
    private final Vector3D gradGamma = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi = (Vector3D)space.makeVector();
    private final Vector3D gradF = (Vector3D)space.makeVector();
    private final Vector3D[] gradEi = new Vector3D[1];
    
    protected NearestImageTransformer nearestImageTransformer;
    protected final Vector rij = (Vector3D)space.makeVector();
    protected final Vector rik = (Vector3D)space.makeVector();
    protected final Vector rkj = (Vector3D)space.makeVector();
    private ParameterSetMEAM p;
    private AtomPair pair;

}
