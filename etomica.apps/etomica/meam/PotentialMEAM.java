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
    	
    	//if (r > 20) return 0;
    	//else {
        
        double sumRhoj0 = 0, sumRhoj1 = 0, sumRhoj2 = 0, sumRhoj3 = 0, 
			sumt1Rhoj0 = 0, sumt2Rhoj0 = 0, sumt3Rhoj0 = 0, 
			sumRhoj1x = 0, sumRhoj1y = 0, sumRhoj1z = 0,  
			sumRhoj2xx = 0, sumRhoj2xy = 0,  sumRhoj2xz = 0,
			sumRhoj2yy = 0, sumRhoj2yz = 0, sumRhoj2zz = 0,
			sumRhoj3xxx = 0, sumRhoj3xxy = 0, sumRhoj3xxz = 0, 
			sumRhoj3xyy = 0, sumRhoj3xyz = 0, sumRhoj3xzz = 0,
			sumRhoj3yyy = 0, sumRhoj3yyz = 0, sumRhoj3yzz = 0, sumRhoj3zzz = 0,
			sumPhi = 0;
		
        
        for(int i = 1; i < atoms.count(); i++) {
        	AtomLeaf atomi = (AtomLeaf) atoms.getAtom(i);
            dr.Ev1Mv2(atomi.coord.position(), atom0.coord.position());
            nearestImageTransformer.nearestImage(dr);
            double r = Math.sqrt(dr.squared());
    		
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
    	double rhoj0 = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0));
    	sumRhoj0 += p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0));
    	double rhoj1 = p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0));
        sumRhoj1 += p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0));
        double rhoj2 = p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0));
        sumRhoj2 += p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0));
    	double rhoj3 = p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0));
    	sumRhoj3 += p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0));
    	
    	sumt1Rhoj0 += rhoj0 * p.t1;
    	sumt2Rhoj0 += rhoj0 * p.t2;
    	sumt3Rhoj0 += rhoj0 * p.t3;
    	
    	unitVector.E(dr);
        unitVector.normalize();
        
        double x = unitVector.x(0);
        double y = unitVector.x(1);
        double z = unitVector.x(2);
       
        sumRhoj1x += rhoj1 * x;
        sumRhoj1y = rhoj1 * y;
        sumRhoj1z = rhoj1 * z;  	
        
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
        double rhoi0Ref = p.Z * rhoj0;
        double rhoi1Ref = p.Z * rhoj1;
        double rhoi2Ref = Math.sqrt(2.0/3.0) * p.Z * rhoj2;
        double rhoi3Ref = p.Z * rhoj3;
        double gammaRef = (p.t1 * (rhoi1Ref/rhoi0Ref) * (rhoi1Ref/rhoi0Ref)) 
			+ (p.t2 * (rhoi2Ref/rhoi0Ref) * (rhoi2Ref/rhoi0Ref)) 
			+ (p.t3 * (rhoi3Ref/rhoi0Ref) * (rhoi3Ref/rhoi0Ref));
        double rhoiRef = (2.0 * rhoi0Ref) / (1.0 + Math.exp(-gammaRef));
        double a = p.alpha * ((r/p.r0) - 1.0);
    	double EuRef = - p.Ec * (1.0 + a) * Math.exp(-a);
    	double FRef = p.A * p.Ec * (rhoiRef / p.rhoScale) * Math.log(rhoiRef/p.rhoScale);
    	sumPhi += (2.0/p.Z) * (EuRef - FRef);
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
        double rhoi = (2.0 * rhoi0) / (1.0 + Math.exp(-gamma));
   
        //The embedding energy of atom i
		double F = p.A * p.Ec * (rhoi/p.rhoScale) * Math.log(rhoi/p.rhoScale);
		
		return F + ((1.0/2.0)*sumPhi);
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
        
        //if (r > 20) return nada;
        //else {
        
        double sumRhoj0 = 0, sumRhoj1 = 0, sumRhoj2 = 0, sumRhoj3 = 0, 
		sumt1Rhoj0 = 0, sumt2Rhoj0 = 0, sumt3Rhoj0 = 0, 
		sumRhoj1x = 0, sumRhoj1y = 0, sumRhoj1z = 0,  
		sumRhoj2xx = 0, sumRhoj2xy = 0,  sumRhoj2xz = 0,
		sumRhoj2yy = 0, sumRhoj2yz = 0, sumRhoj2zz = 0,
		sumRhoj3xxx = 0, sumRhoj3xxy = 0, sumRhoj3xxz = 0, 
		sumRhoj3xyy = 0, sumRhoj3xyz = 0, sumRhoj3xzz = 0,
		sumRhoj3yyy = 0, sumRhoj3yyz = 0, sumRhoj3yzz = 0, sumRhoj3zzz = 0,
		sumPhi = 0;
        
        for(int i = 1; i < atoms.count(); i++) {
        	
        	AtomLeaf atomi = (AtomLeaf) atoms.getAtom(i);
            dr.Ev1Mv2(atomi.coord.position(), atom0.coord.position());
            nearestImageTransformer.nearestImage(dr);
            double r = Math.sqrt(dr.squared());
        	
        	double rhoj0 = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0));
        	sumRhoj0 += p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1.0));
        	double rhoj1 = p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0));
            sumRhoj1 += p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1.0));
            double rhoj2 = p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0));
            sumRhoj2 += p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1.0));
        	double rhoj3 = p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0));
        	sumRhoj3 += p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1.0));
        	
        	sumt1Rhoj0 += rhoj0 * p.t1;
        	sumt2Rhoj0 += rhoj0 * p.t2;
        	sumt3Rhoj0 += rhoj0 * p.t3;
        	
        	unitVector.E(dr);
            unitVector.normalize();
            
            double x = unitVector.x(0);
            double y = unitVector.x(1);
            double z = unitVector.x(2);
           
            sumRhoj1x += rhoj1 * x;
            sumRhoj1y = rhoj1 * y;
            sumRhoj1z = rhoj1 * z;  	
            
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
        
    	double rhoi0Ref = p.Z * rhoj0;
        double rhoi1Ref = p.Z * rhoj1;
        double rhoi2Ref = Math.sqrt(2.0/3.0) * p.Z * rhoj2;
        double rhoi3Ref = p.Z * rhoj3;
        double gammaRef = (p.t1 * (rhoi1Ref/rhoi0Ref) * (rhoi1Ref/rhoi0Ref)) 
 			+ (p.t2 * (rhoi2Ref/rhoi0Ref) * (rhoi2Ref/rhoi0Ref)) 
 			+ (p.t3 * (rhoi3Ref/rhoi0Ref) * (rhoi3Ref/rhoi0Ref));
         double rhoiRef = (2.0 * rhoi0Ref) / (1.0 + Math.exp(-gammaRef));
         
    	
    	vector100.setX(0,1);
    	vector010.setX(1,1);
    	vector001.setX(2,1);
    	
    	gradx.Ea1Tv1(-x/(r*r),dr);
    	gradx.PEa1Tv1(1/r, vector100);
    	
    	grady.Ea1Tv1(-y/(r*r),dr);
    	grady.PEa1Tv1(1/r, vector010);
    	
    	gradz.Ea1Tv1(-z/(r*r),dr);
    	gradz.PEa1Tv1(1/r, vector001);
    	
    	gradr.Ea1Tv1(1/r,dr);
    	
    	//Gradient of rhoj0
    	gradRhoj0.Ea1Tv1(-rhoj0*p.beta0/(p.r0), gradr);
    	sumGradRhoj0.PE(gradRhoj0);
    	
    	//Gradient of rhoj1
    	gradRhoj1.Ea1Tv1(-rhoj1*p.beta1/(p.r0), gradr);
    	sumGradRhoj1.PE(gradRhoj1);
    	
    	//Gradient of rhoj2
    	gradRhoj2.Ea1Tv1(-rhoj2*p.beta2/(p.r0), gradr);
    	sumGradRhoj2.PE(gradRhoj2);
    	
    	//Gradient of rhoj3
    	gradRhoj3.Ea1Tv1(-rhoj3*p.beta3/(p.r0), gradr);
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
    	gradRhoj3xxx.Ea1Tv1(6.0*rhoj3*x*x, gradx);
    	gradRhoj3xxx.PEa1Tv1(x*x*x, gradRhoj3);
    	sumGradRhoj3xxx.PE(gradRhoj3xxx);
    	
    	//Gradient of rhoj3xxy
    	gradRhoj3xxy.Ea1Tv1(2.0*rhoj3*x*x, grady);
    	gradRhoj3xxy.PEa1Tv1(4.0*rhoj3*x*y, gradx);
    	gradRhoj3xxy.PEa1Tv1(x*x*y, gradRhoj3);
    	sumGradRhoj3xxy.PE(gradRhoj3xxy);
    	
    	//Gradient of rhoj3xxz
    	gradRhoj3xxz.Ea1Tv1(2.0*rhoj3*x*x, gradz);
    	gradRhoj3xxz.PEa1Tv1(4.0*rhoj3*x*z, gradx);
    	gradRhoj3xxz.PEa1Tv1(x*x*z, gradRhoj3);
    	sumGradRhoj3xxz.PE(gradRhoj3xxz);
    	
    	//Gradient of rhoj3xyy
    	gradRhoj3xyy.Ea1Tv1(4.0*rhoj3*x*y, grady);
    	gradRhoj3xyy.PEa1Tv1(2.0*rhoj3*y*y, gradx);
    	gradRhoj3xyy.PEa1Tv1(x*y*y, gradRhoj3);
    	sumGradRhoj3xyy.PE(gradRhoj3xyy);
    	
    	//Gradient of rhoj3xyz
    	gradRhoj3xyz.Ea1Tv1(2.0*rhoj3*x*y, gradz);
    	gradRhoj3xyz.PEa1Tv1(2.0*rhoj3*x*z, grady);
    	gradRhoj3xyz.PEa1Tv1(2.0*rhoj3*y*z, gradx);
    	gradRhoj3xyz.PEa1Tv1(x*y*z, gradRhoj3);
    	sumGradRhoj3xyz.PE(gradRhoj3xyz);
    	
    	//Gradient of rhoj3xzz
    	gradRhoj3xzz.Ea1Tv1(4.0*rhoj3*x*z, gradz);
    	gradRhoj3xzz.PEa1Tv1(2.0*rhoj3*z*z, gradx);
    	gradRhoj3xzz.PEa1Tv1(x*z*z, gradRhoj3);
    	sumGradRhoj3xzz.PE(gradRhoj3xzz);
    	
    	//Gradient of rhoj3yyy
    	gradRhoj3yyy.Ea1Tv1(6.0*rhoj3*y*y, grady);
    	gradRhoj3yyy.PEa1Tv1(y*y*y, gradRhoj3);
    	sumGradRhoj3yyy.PE(gradRhoj3yyy);
    	
    	//Gradient of rhoj3yyz
    	gradRhoj3yyz.Ea1Tv1(2.0*rhoj3*y*y, gradz);
    	gradRhoj3yyz.PEa1Tv1(4.0*rhoj3*y*z, grady);
    	gradRhoj3yyz.PEa1Tv1(y*y*z, gradRhoj3);
    	sumGradRhoj3yyz.PE(gradRhoj3yyz);
    	
    	//Gradient of rhoj3yzz
    	gradRhoj3yzz.Ea1Tv1(4.0*rhoj3*y*z, gradz);
    	gradRhoj3yzz.PEa1Tv1(2.0*rhoj3*z*z, grady);
    	gradRhoj3yzz.PEa1Tv1(y*z*z, gradRhoj3);
    	sumGradRhoj3yzz.PE(gradRhoj3yzz);
    	
    	//Gradient of rhoj3zzz
    	gradRhoj3zzz.Ea1Tv1(6.0*rhoj3*z*z, gradz);
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
    	
    	//Calculations to determine gradient of phi
    	
    	gradRhoi0Ref.Ea1Tv1(p.Z, gradRhoj0);
    	
    	gradRhoi1Ref.Ea1Tv1(x, gradRhoj1x);
    	gradRhoi1Ref.PEa1Tv1(y, gradRhoj1y);
    	gradRhoi1Ref.PEa1Tv1(z, gradRhoj1z);
    	gradRhoi1Ref.TE(p.Z*p.Z*rhoj1/rhoi1Ref);
    	
    	gradRhoi2Ref.Ea1Tv1(x*x, gradRhoj2xx);
    	gradRhoi2Ref.PEa1Tv1(2.0*x*y, gradRhoj2xy);
    	gradRhoi2Ref.PEa1Tv1(2.0*x*z, gradRhoj2xz);
    	gradRhoi2Ref.PEa1Tv1(y*y, gradRhoj2yy);
    	gradRhoi2Ref.PEa1Tv1(2.0*y*z, gradRhoj2yz);
    	gradRhoi2Ref.PEa1Tv1(z*z, gradRhoj2zz);
    	gradRhoi2Ref.PEa1Tv1(-(1.0/3.0), gradRhoj2);
    	gradRhoi2Ref.TE(p.Z*p.Z*rhoj2/rhoi2Ref);
    	
    	gradRhoi3Ref.Ea1Tv1(x*x*x, gradRhoj3xxx);
    	gradRhoi3Ref.PEa1Tv1(3.0*x*x*y, gradRhoj3xxy);
    	gradRhoi3Ref.PEa1Tv1(3.0*x*x*z, gradRhoj3xxz);
    	gradRhoi3Ref.PEa1Tv1(3.0*x*y*y, gradRhoj3xyy);
    	gradRhoi3Ref.PEa1Tv1(6.0*x*y*z, gradRhoj3xyz);
    	gradRhoi3Ref.PEa1Tv1(3.0*x*z*z, gradRhoj3xzz);
    	gradRhoi3Ref.PEa1Tv1(y*y*y, gradRhoj3yyy);
    	gradRhoi3Ref.PEa1Tv1(3.0*y*y*z, gradRhoj3yyz);
    	gradRhoi3Ref.PEa1Tv1(3.0*y*z*z, gradRhoj3yzz);
    	gradRhoi3Ref.PEa1Tv1(z*z*z, gradRhoj3zzz);
    	gradRhoi3Ref.TE(p.Z*p.Z*rhoj3/rhoi3Ref);
    	
    	gradGammaRef.Ea1Tv1( (-1.0/rhoi0Ref) * (p.t1*rhoi1Ref*rhoi1Ref +
    			p.t2*rhoi2Ref*rhoi2Ref + p.t3*rhoi3Ref*rhoi3Ref), gradRhoi0Ref);
    	gradGammaRef.PEa1Tv1(p.t1*rhoi1Ref, gradRhoi1Ref);
    	gradGammaRef.PEa1Tv1(p.t2*rhoi2Ref, gradRhoi2Ref);
    	gradGammaRef.PEa1Tv1(p.t3*rhoi3Ref, gradRhoi3Ref);
    	gradGammaRef.TE(2.0/(rhoi0Ref*rhoi0Ref));
    	
    	gradRhoiRef.Ea1Tv1(rhoi0Ref*Math.exp(-gammaRef)
    			/((1.0+Math.exp(-gammaRef))*(1.0+Math.exp(-gammaRef))), gradGammaRef);
    	gradRhoiRef.PEa1Tv1(1.0/(1.0+Math.exp(-gammaRef)), gradRhoi0Ref);
    	gradRhoiRef.TE(2.0);
    	
    	gradFRef.Ea1Tv1( (p.A*p.Ec / p.rhoScale)*
    			(1.0 + Math.log(rhoiRef) - Math.log(p.rhoScale)), gradRhoiRef);
    	
    	gradERef.Ea1Tv1( (p.Ec*p.alpha*p.alpha/p.r0) * ((r/p.r0) - 1.0)
    					*(Math.exp(-p.alpha*((r/p.r0)-1.0))), gradr);
    	
    	gradPhi.E(gradERef);
    	gradPhi.ME(gradFRef);
    	gradPhi.TE(2.0/p.Z);
    	sumGradPhi.PE(gradPhi);
    	
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
        double rhoi = (2.0 * rhoi0) / (1.0 + Math.exp(-gamma));
        
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
		
		gradRhoi.Ea1Tv1(rhoi0*Math.exp(-gamma)/(1.0 + Math.exp(-gamma)), gradGamma);
		gradRhoi.PE(gradRhoi0);
		gradRhoi.TE(2.0/(1.0+Math.exp(-gamma)));
	
		gradF.Ea1Tv1( (p.A*p.Ec/p.rhoScale)*
				(1.0 + Math.log(rhoi)- Math.log(p.rhoScale)), gradRhoi);
		
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
    private final Vector3D gradRhoi0Ref = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi1Ref = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi2Ref = (Vector3D)space.makeVector();
    private final Vector3D gradRhoi3Ref = (Vector3D)space.makeVector();
    private final Vector3D gradGammaRef = (Vector3D)space.makeVector();
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
    protected final Vector dr = (Vector3D)space.makeVector();
    private ParameterSetMEAM p;
    private AtomPair pair;

}
