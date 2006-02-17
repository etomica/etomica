package etomica.meam;
import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.phase.Phase;
import etomica.phase.PhaseAgentManager;
import etomica.phase.PhaseAgentSourceAtomManager;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.simulation.Simulation;
import etomica.space.CoordinatePair;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.Dimension;
import etomica.units.Energy;

/**
 * This class calculates the pair-wise terms required for the modified embedded-
 * atom method potential.  The values calculated that are relevant to each atom
 * are stored and summed in bins, which may be accessed later when calculating 
 * each atom's energy.
 */

public final class MEAMP2 extends Potential2 implements Potential2Soft, AtomAgentManager.AgentSource {

	public MEAMP2(Simulation sim, ParameterSetMEAM p) {
		super(sim.space);
        this.p = p;
        coordinatePair = new CoordinatePair(sim.space); 
        unitVector = (Vector3D)space.makeVector();
        vector100 = (Vector3D)space.makeVector();
        vector010 = (Vector3D)space.makeVector();
        vector001 = (Vector3D)space.makeVector();
        gradx = (Vector3D)space.makeVector();
        grady = (Vector3D)space.makeVector();
        gradz = (Vector3D)space.makeVector();
        gradr = (Vector3D)space.makeVector();
        gradRhoj0 = (Vector3D)space.makeVector();
        gradRhoj1 = (Vector3D)space.makeVector();
        gradRhoj2 = (Vector3D)space.makeVector();
        gradRhoj3 = (Vector3D)space.makeVector();
        gradRhoj1x = (Vector3D)space.makeVector();
        gradRhoj1y = (Vector3D)space.makeVector();
        gradRhoj1z = (Vector3D)space.makeVector();
        gradRhoj2xx = (Vector3D)space.makeVector();
        gradRhoj2xy = (Vector3D)space.makeVector();
        gradRhoj2xz = (Vector3D)space.makeVector();
        gradRhoj2yy = (Vector3D)space.makeVector();
        gradRhoj2yz = (Vector3D)space.makeVector();
        gradRhoj2zz = (Vector3D)space.makeVector();
        gradRhoj3xxx = (Vector3D)space.makeVector();
        gradRhoj3xxy = (Vector3D)space.makeVector();
        gradRhoj3xxz = (Vector3D)space.makeVector();
        gradRhoj3xyy = (Vector3D)space.makeVector();
        gradRhoj3xyz = (Vector3D)space.makeVector();
        gradRhoj3xzz = (Vector3D)space.makeVector();
        gradRhoj3yyy = (Vector3D)space.makeVector();
        gradRhoj3yyz = (Vector3D)space.makeVector();
        gradRhoj3yzz = (Vector3D)space.makeVector();
        gradRhoj3zzz = (Vector3D)space.makeVector();
        t1GradRhoj0 = (Vector3D)space.makeVector();
        t2GradRhoj0 = (Vector3D)space.makeVector();
        t3GradRhoj0 = (Vector3D)space.makeVector();
        gradRhoi0Ref = (Vector3D)space.makeVector();
        gradRhoi1Ref = (Vector3D)space.makeVector();
        gradRhoi2Ref = (Vector3D)space.makeVector();
        gradRhoj213 = (Vector3D)space.makeVector();
        gradRhoi3Ref = (Vector3D)space.makeVector();
        gradGammaRef = (Vector3D)space.makeVector();
        gradRhoiRef = (Vector3D)space.makeVector();
        gradERef = (Vector3D)space.makeVector();
        gradFRef = (Vector3D)space.makeVector();
        gradPhi = (Vector3D)space.makeVector();
        nada = (Vector3D)space.makeVector();
        
        
        
        phaseAgentManager = new PhaseAgentManager(new PhaseAgentSourceAtomManager(this),sim.speciesRoot);
    }
	
    public void setPhase(Phase phase){
    	coordinatePair.setNearestImageTransformer(phase.getBoundary());
        agentManager = (AtomAgentManager[])phaseAgentManager.getAgents();
        agents = (Wrapper[])agentManager[phase.getIndex()].getAgents();
    }

    /**
     * The following field calculates the values of the pair-wise
     * terms required for calculation of the embedding energy of both atoms in 
     * the pair, and it also calculates the value of phi, the repulsive energy
     * between each pair of atoms, another pair-wise term.  
     * 
     * The values required for each atom are stored in separate bins for
     * each atom in an array.  When another pair is examined that includes an
     * atom already examined in a previous pairing, the new pair-wise terms 
     * calculated for the new pairing will be added to the previous values.  In
     * this way, the summations over all neighbor's j are performed.   
     * 
     * This field functions only for unmixed pairs of atoms.  The parameters, 
     * which are dependent upon the type of elements involved in the pair, 
     * are provided in the class ParameterSetEAM.
     */
    public double energy(AtomSet pair) {
        coordinatePair.reset((AtomPair)pair);
    	double r = Math.sqrt(coordinatePair.r2());
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
    	double rhoj0 = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1));
    	double rhoj1 = p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1));
    	double rhoj2 = p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1));
    	double rhoj3 = p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1));
    	//
    	//We need to store the information required by the many-body potential 
    	//so that it may be accessed when MEAMPMany is employed.  Information 
    	//relevant to each atom in the pair of atoms (i and j, 0 and 1), 
    	//currently being examined is sent to an array.  The NEAMPMany class then
    	//only needs to access one atom's information in the array to determine
    	//the embedding energy of the atom.
    	//
    	//The information being stored is that which is required before the 
    	//summations over neighbors j in the expressions for partial electron 
    	//densities of atom i.
    	//
    	double[] sumsAtom0 = agents[((AtomPair)pair).atom0.getGlobalIndex()].sums;
    	double[] sumsAtom1 = agents[((AtomPair)pair).atom1.getGlobalIndex()].sums;
    	
    	double t1Rhoj0 = rhoj0 * p.t1;
    	double t2Rhoj0 = rhoj0 * p.t2;
    	double t3Rhoj0 = rhoj0 * p.t3;
    	
    	sumsAtom0[Wrapper.rhoj0Bin] += rhoj0;
    	sumsAtom1[Wrapper.rhoj0Bin] += rhoj0;
    	sumsAtom0[Wrapper.t1Rhoj0Bin] += t1Rhoj0;
    	sumsAtom1[Wrapper.t1Rhoj0Bin] += t1Rhoj0;
    	sumsAtom0[Wrapper.t2Rhoj0Bin] += t2Rhoj0;
    	sumsAtom1[Wrapper.t2Rhoj0Bin] += t2Rhoj0;
    	sumsAtom0[Wrapper.t3Rhoj0Bin] += t3Rhoj0;
    	sumsAtom1[Wrapper.t3Rhoj0Bin] += t3Rhoj0;
    	
        unitVector.E(coordinatePair.dr());
        unitVector.normalize();
        
        double x = unitVector.x(0);
        double y = unitVector.x(1);
        double z = unitVector.x(2);
        
        double rhoj1x = rhoj1 * x;
        double rhoj1y = rhoj1 * y;
        double rhoj1z = rhoj1 * z;
        
        
        sumsAtom0[Wrapper.rhoj1Bin] += rhoj1;
        sumsAtom1[Wrapper.rhoj1Bin] += rhoj1;
        sumsAtom0[Wrapper.rhoj1xBin] += rhoj1x;
        sumsAtom1[Wrapper.rhoj1xBin] += rhoj1x;
        sumsAtom0[Wrapper.rhoj1yBin] += rhoj1y;
        sumsAtom1[Wrapper.rhoj1yBin] += rhoj1y;
        sumsAtom0[Wrapper.rhoj1zBin] += rhoj1z;
        sumsAtom1[Wrapper.rhoj1zBin] += rhoj1z;
        
        double rhoj2xx = rhoj2 * x * x;
        double rhoj2xy = rhoj2 * x * y;
        double rhoj2xz = rhoj2 * x * z;
        double rhoj2yy = rhoj2 * y * y;
        double rhoj2yz = rhoj2 * y * z;
        double rhoj2zz = rhoj2 * z * z;
        
        
        sumsAtom0[Wrapper.rhoj2Bin] += rhoj2;
        sumsAtom1[Wrapper.rhoj2Bin] += rhoj2;
        sumsAtom0[Wrapper.rhoj2xxBin] += rhoj2xx;
        sumsAtom1[Wrapper.rhoj2xxBin] += rhoj2xx;
        sumsAtom0[Wrapper.rhoj2xyBin] += rhoj2xy;
        sumsAtom1[Wrapper.rhoj2xyBin] += rhoj2xy;
        sumsAtom0[Wrapper.rhoj2xzBin] += rhoj2xz;
        sumsAtom1[Wrapper.rhoj2xzBin] += rhoj2xz;
        sumsAtom0[Wrapper.rhoj2yyBin] += rhoj2yy;
        sumsAtom1[Wrapper.rhoj2yyBin] += rhoj2yy;
        sumsAtom0[Wrapper.rhoj2yzBin] += rhoj2yz;
        sumsAtom1[Wrapper.rhoj2yzBin] += rhoj2yz;
        sumsAtom0[Wrapper.rhoj2zzBin] += rhoj2zz;
        sumsAtom1[Wrapper.rhoj2zzBin] += rhoj2zz;
        
        double rhoj3xxx = rhoj3 * x * x * x;
        double rhoj3xxy = rhoj3 * x * x * y;
        double rhoj3xxz = rhoj3 * x * x * z;
        double rhoj3xyy = rhoj3 * x * y * y;
        double rhoj3xyz = rhoj3 * x * y * z;
        double rhoj3xzz = rhoj3 * x * z * z;
        double rhoj3yyy = rhoj3 * y * y * y;
        double rhoj3yyz = rhoj3 * y * y * z;
        double rhoj3yzz = rhoj3 * y * z * z;
        double rhoj3zzz = rhoj3 * z * z * z;
        
        sumsAtom0[Wrapper.rhoj3Bin] += rhoj3;
        sumsAtom1[Wrapper.rhoj3Bin] += rhoj3;
        sumsAtom0[Wrapper.rhoj3xxxBin] += rhoj3xxx;
        sumsAtom1[Wrapper.rhoj3xxxBin] += rhoj3xxx;
        sumsAtom0[Wrapper.rhoj3xxyBin] += rhoj3xxy;
        sumsAtom1[Wrapper.rhoj3xxyBin] += rhoj3xxy;
        sumsAtom0[Wrapper.rhoj3xxzBin] += rhoj3xxz;
        sumsAtom1[Wrapper.rhoj3xxzBin] += rhoj3xxz;
        sumsAtom0[Wrapper.rhoj3xyyBin] += rhoj3xyy;
        sumsAtom1[Wrapper.rhoj3xyyBin] += rhoj3xyy;
        sumsAtom0[Wrapper.rhoj3xyzBin] += rhoj3xyz;
        sumsAtom1[Wrapper.rhoj3xyzBin] += rhoj3xyz;
        sumsAtom0[Wrapper.rhoj3xzzBin] += rhoj3xzz;
        sumsAtom1[Wrapper.rhoj3xzzBin] += rhoj3xzz;
        sumsAtom0[Wrapper.rhoj3yyyBin] += rhoj3yyy;
        sumsAtom1[Wrapper.rhoj3yyyBin] += rhoj3yyy;
        sumsAtom0[Wrapper.rhoj3yyzBin] += rhoj3yyz;
        sumsAtom1[Wrapper.rhoj3yyzBin] += rhoj3yyz;
        sumsAtom0[Wrapper.rhoj3yzzBin] += rhoj3yzz;
        sumsAtom1[Wrapper.rhoj3yzzBin] += rhoj3yzz;
        sumsAtom0[Wrapper.rhoj3zzzBin] += rhoj3zzz;
        sumsAtom1[Wrapper.rhoj3zzzBin] += rhoj3zzz;

        //
        //Now, we continue to calculate the pair potential, phi.
        //
        double rhoi0Ref = p.Z * rhoj0;
        double rhoi1Ref = p.Z * rhoj1;
        double rhoi2Ref = Math.sqrt(2/3) * p.Z * rhoj2;
        double rhoi3Ref = p.Z * rhoj3;
        double gammaRef = (p.t1 * (rhoi1Ref/rhoi0Ref) * (rhoi1Ref/rhoi0Ref)) 
			+ (p.t2 * (rhoi2Ref/rhoi0Ref) * (rhoi2Ref/rhoi0Ref)) 
			+ (p.t3 * (rhoi3Ref/rhoi0Ref) * (rhoi3Ref/rhoi0Ref));
        double rhoiRef = (2 * rhoi0Ref) / (1 + Math.exp(-gammaRef));
        double a = p.alpha * ((r/p.r0) - 1);
    	double EuRef = - p.Ec * (1 + a) * Math.exp(- a);
    	double FRef = p.A * p.Ec * (rhoiRef / p.rhoScale) 
			* (Math.log(rhoiRef / p.rhoScale));
    	double phi = (2/p.Z) * (EuRef - FRef);
    	
    	//We also need to store phi for each atom, so that it may be summed 
    	//as all of the pairs are evaluated.
    	
    	sumsAtom0[Wrapper.phiBin] += phi;
    	sumsAtom1[Wrapper.phiBin] += phi;
    	
    	//MEAMP2 should contribute nothing directly to the potential.
    	return 0;
    }
    
    /**
     * Virial of the pair as given by the du(double) method
     */
    public double virial(AtomSet pair) {
    	coordinatePair.reset((AtomPair)pair);
        return 0;
    }
    
    /**
     * Hypervirial of the pair as given by the du(double) and d2u(double) methods
     */
    public double hyperVirial(AtomSet pair) {
    	coordinatePair.reset((AtomPair)pair);
        double r2 = coordinatePair.r2();
        return 0;
    }
    
    /**
     * Gradient of the pair potential as given by the du(double) method.
     */
    public Vector gradient(AtomSet pair) {
    	coordinatePair.reset((AtomPair)pair);
        double r2 = coordinatePair.r2();
        double r = Math.sqrt(r2);
        this.pair = (AtomPair)pair;
        
        double rhoj0 = p.rhoScale * Math.exp(-p.beta0 * ((r/p.r0) - 1));
    	double rhoj1 = p.rhoScale * Math.exp(-p.beta1 * ((r/p.r0) - 1));
    	double rhoj2 = p.rhoScale * Math.exp(-p.beta2 * ((r/p.r0) - 1));
    	double rhoj3 = p.rhoScale * Math.exp(-p.beta3 * ((r/p.r0) - 1));
    	double rhoi0Ref = p.Z * rhoj0;
        double rhoi1Ref = p.Z * rhoj1;
        double rhoi2Ref = Math.sqrt(2/3) * p.Z * rhoj2;
        double rhoi3Ref = p.Z * rhoj3;
        double gammaRef = (p.t1 * (rhoi1Ref/rhoi0Ref) * (rhoi1Ref/rhoi0Ref)) 
 			+ (p.t2 * (rhoi2Ref/rhoi0Ref) * (rhoi2Ref/rhoi0Ref)) 
 			+ (p.t3 * (rhoi3Ref/rhoi0Ref) * (rhoi3Ref/rhoi0Ref));
         double rhoiRef = (2 * rhoi0Ref) / (1 + Math.exp(-gammaRef));
    	
    	unitVector.E(coordinatePair.dr());
        unitVector.normalize();
        
        double x = unitVector.x(0);
        double y = unitVector.x(1);
        double z = unitVector.x(2);
        
//      Calculations of terms required for gradient of potential energy
    	
    	Vector3D[] gradientSumsAtom0 = agents[((AtomPair)pair).atom0.getGlobalIndex()].gradientSums;
    	Vector3D[] gradientSumsAtom1 = agents[((AtomPair)pair).atom1.getGlobalIndex()].gradientSums;
    	
    	vector100.setX(0,1);
    	vector010.setX(1,1);
    	vector001.setX(2,1);
    	
    	gradx.Ea1Tv1(-x/(r*r),coordinatePair.dr());
    	gradx.PEa1Tv1(1/r, vector100);
    	grady.Ea1Tv1(-y/(r*r),coordinatePair.dr());
    	grady.PEa1Tv1(1/r, vector010);
    	gradz.Ea1Tv1(-z/(r*r),coordinatePair.dr());
    	gradz.PEa1Tv1(1/r, vector001);
    	
    	gradr.Ea1Tv1(1/r,coordinatePair.dr());
    	
    	//Gradient of rhoj0
    	gradRhoj0.Ea1Tv1(-rhoj0*p.beta0/(p.r0), gradr);
    	
    	gradientSumsAtom0[Wrapper.rhoj0Bin].PE(gradRhoj0);
    	gradientSumsAtom1[Wrapper.rhoj0Bin].PE(gradRhoj0);
    	
    	//Gradient of rhoj1
    	gradRhoj1.Ea1Tv1(-rhoj1*p.beta1/(p.r0), gradr);
    	
    	gradientSumsAtom0[Wrapper.rhoj1Bin].PE(gradRhoj1);
    	gradientSumsAtom1[Wrapper.rhoj1Bin].PE(gradRhoj1);
    	
    	//Gradient of rhoj2
    	gradRhoj2.Ea1Tv1(-rhoj2*p.beta2/(p.r0), gradr);
    	
    	gradientSumsAtom0[Wrapper.rhoj2Bin].PE(gradRhoj2);
    	gradientSumsAtom1[Wrapper.rhoj2Bin].PE(gradRhoj2);
    	
    	//Gradient of rhoj3
    	gradRhoj3.Ea1Tv1(-rhoj3*p.beta3/(p.r0), gradr);
    	
    	gradientSumsAtom0[Wrapper.rhoj3Bin].PE(gradRhoj3);
    	gradientSumsAtom1[Wrapper.rhoj3Bin].PE(gradRhoj3);
    	
    	//Gradient of rhoj1x
    	gradRhoj1x.Ea1Tv1(rhoj1, gradx);
    	gradRhoj1x.PEa1Tv1(x, gradRhoj1);
    	
    	gradientSumsAtom0[Wrapper.rhoj1xBin].PE(gradRhoj1x);
    	gradientSumsAtom1[Wrapper.rhoj1xBin].PE(gradRhoj1x);
    	
    	//Gradient of rhoj1y
    	gradRhoj1y.Ea1Tv1(rhoj1, grady);
    	gradRhoj1y.PEa1Tv1(y, gradRhoj1);
    	
    	gradientSumsAtom0[Wrapper.rhoj1yBin].PE(gradRhoj1y);
    	gradientSumsAtom1[Wrapper.rhoj1yBin].PE(gradRhoj1y);
    	
    	//Gradient of rhoj1z
    	gradRhoj1z.Ea1Tv1(rhoj1, gradz);
    	gradRhoj1z.PEa1Tv1(z, gradRhoj1);
    	
    	gradientSumsAtom0[Wrapper.rhoj1zBin].PE(gradRhoj1z);
    	gradientSumsAtom1[Wrapper.rhoj1zBin].PE(gradRhoj1z);
    	
    	//Gradient of rhoj2xx
    	gradRhoj2xx.Ea1Tv1(2*rhoj2*x, gradx);
    	gradRhoj2xx.PEa1Tv1(x*x, gradRhoj2);
    	
    	gradientSumsAtom0[Wrapper.rhoj2xxBin].PE(gradRhoj2xx);
    	gradientSumsAtom1[Wrapper.rhoj2xxBin].PE(gradRhoj2xx);
    	
    	//Gradient of rhoj2xy
    	gradRhoj2xy.Ea1Tv1(rhoj2*x, grady);
    	gradRhoj2xy.Ea1Tv1(rhoj2*y, gradx);
    	gradRhoj2xy.PEa1Tv1(x*y, gradRhoj2);
    	
    	gradientSumsAtom0[Wrapper.rhoj2xyBin].PE(gradRhoj2xy);
    	gradientSumsAtom1[Wrapper.rhoj2xyBin].PE(gradRhoj2xy);
    	
    	//Gradient of rhoj2xz
    	gradRhoj2xz.Ea1Tv1(rhoj2*x, gradz);
    	gradRhoj2xz.Ea1Tv1(rhoj2*z, gradx);
    	gradRhoj2xz.PEa1Tv1(x*z, gradRhoj2);
    	
    	gradientSumsAtom0[Wrapper.rhoj2xzBin].PE(gradRhoj2xz);
    	gradientSumsAtom1[Wrapper.rhoj2xzBin].PE(gradRhoj2xz);
    	
    	//Gradient of rhoj2yy
    	gradRhoj2yy.Ea1Tv1(2*rhoj2*y, grady);
    	gradRhoj2yy.PEa1Tv1(y*y, gradRhoj2);
    	
    	gradientSumsAtom0[Wrapper.rhoj2yyBin].PE(gradRhoj2yy);
    	gradientSumsAtom1[Wrapper.rhoj2yyBin].PE(gradRhoj2yy);
    	
    	//Gradient of rhoj2yz
    	gradRhoj2yz.Ea1Tv1(rhoj2*y, gradz);
    	gradRhoj2yz.Ea1Tv1(rhoj2*z, grady);
    	gradRhoj2yz.PEa1Tv1(y*z, gradRhoj2);
    	
    	gradientSumsAtom0[Wrapper.rhoj2yzBin].PE(gradRhoj2yz);
    	gradientSumsAtom1[Wrapper.rhoj2yzBin].PE(gradRhoj2yz);
    	
    	//Gradient of rhoj2zz
    	gradRhoj2zz.Ea1Tv1(2*rhoj2*z, gradz);
    	gradRhoj2zz.PEa1Tv1(z*z, gradRhoj2);
    	
    	gradientSumsAtom0[Wrapper.rhoj2zzBin].PE(gradRhoj2zz);
    	gradientSumsAtom1[Wrapper.rhoj2zzBin].PE(gradRhoj2zz);
    	
    	//Gradient of rhoj3xxx
    	gradRhoj3xxx.Ea1Tv1(6*rhoj3*x*x, gradx);
    	gradRhoj3xxx.PEa1Tv1(x*x*x, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3xxxBin].PE(gradRhoj3xxx);
    	gradientSumsAtom1[Wrapper.rhoj3xxxBin].PE(gradRhoj3xxx);
    	
    	//Gradient of rhoj3xxy
    	gradRhoj3xxy.Ea1Tv1(2*rhoj3*x*x, grady);
    	gradRhoj3xxy.PEa1Tv1(4*rhoj3*x*y, gradx);
    	gradRhoj3xxy.PEa1Tv1(x*x*y, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3xxyBin].PE(gradRhoj3xxy);
    	gradientSumsAtom1[Wrapper.rhoj3xxyBin].PE(gradRhoj3xxy);
    	
    	//Gradient of rhoj3xxz
    	gradRhoj3xxz.Ea1Tv1(2*rhoj3*x*x, gradz);
    	gradRhoj3xxz.PEa1Tv1(4*rhoj3*x*z, gradx);
    	gradRhoj3xxz.PEa1Tv1(x*x*z, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3xxzBin].PE(gradRhoj3xxz);
    	gradientSumsAtom1[Wrapper.rhoj3xxzBin].PE(gradRhoj3xxz);
    	
    	//Gradient of rhoj3xyy
    	gradRhoj3xyy.Ea1Tv1(4*rhoj3*x*y, grady);
    	gradRhoj3xyy.PEa1Tv1(2*rhoj3*y*y, gradx);
    	gradRhoj3xyy.PEa1Tv1(x*y*y, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3xyyBin].PE(gradRhoj3xyy);
    	gradientSumsAtom1[Wrapper.rhoj3xyyBin].PE(gradRhoj3xyy);
    	
    	//Gradient of rhoj3xyz
    	gradRhoj3xyz.Ea1Tv1(2*rhoj3*x*y, gradz);
    	gradRhoj3xyz.PEa1Tv1(2*rhoj3*x*z, grady);
    	gradRhoj3xyz.PEa1Tv1(2*rhoj3*y*z, gradx);
    	gradRhoj3xyz.PEa1Tv1(x*y*z, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3xyzBin].PE(gradRhoj3xyz);
    	gradientSumsAtom1[Wrapper.rhoj3xyzBin].PE(gradRhoj3xyz);
    	
    	//Gradient of rhoj3xzz
    	gradRhoj3xzz.Ea1Tv1(4*rhoj3*x*z, gradz);
    	gradRhoj3xzz.PEa1Tv1(2*rhoj3*z*z, gradx);
    	gradRhoj3xzz.PEa1Tv1(x*z*z, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3xzzBin].PE(gradRhoj3xzz);
    	gradientSumsAtom1[Wrapper.rhoj3xzzBin].PE(gradRhoj3xzz);
    	
    	//Gradient of rhoj3yyy
    	gradRhoj3yyy.Ea1Tv1(6*rhoj3*y*y, grady);
    	gradRhoj3yyy.PEa1Tv1(y*y*y, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3yyyBin].PE(gradRhoj3yyy);
    	gradientSumsAtom1[Wrapper.rhoj3yyyBin].PE(gradRhoj3yyy);
    	
    	//Gradient of rhoj3yyz
    	gradRhoj3yyz.Ea1Tv1(2*rhoj3*y*y, gradz);
    	gradRhoj3yyz.PEa1Tv1(4*rhoj3*y*z, grady);
    	gradRhoj3xxz.PEa1Tv1(y*y*z, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3yyzBin].PE(gradRhoj3yyz);
    	gradientSumsAtom1[Wrapper.rhoj3yyzBin].PE(gradRhoj3yyz);
    	
    	//Gradient of rhoj3yzz
    	gradRhoj3yzz.Ea1Tv1(4*rhoj3*y*z, gradz);
    	gradRhoj3yzz.PEa1Tv1(2*rhoj3*z*z, grady);
    	gradRhoj3yzz.PEa1Tv1(y*z*z, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3yzzBin].PE(gradRhoj3yzz);
    	gradientSumsAtom1[Wrapper.rhoj3yzzBin].PE(gradRhoj3yzz);
    	
    	//Gradient of rhoj3zzz
    	gradRhoj3zzz.Ea1Tv1(6*rhoj3*z*z, gradz);
    	gradRhoj3zzz.PEa1Tv1(z*z*z, gradRhoj3);
    	
    	gradientSumsAtom0[Wrapper.rhoj3zzzBin].PE(gradRhoj3zzz);
    	gradientSumsAtom1[Wrapper.rhoj3zzzBin].PE(gradRhoj3zzz);
    	
    	//t1GradRhoj0
    	t1GradRhoj0.Ea1Tv1(p.t1, gradRhoj0);
    	
    	gradientSumsAtom0[Wrapper.t1Rhoj0Bin].PE(t1GradRhoj0);
    	gradientSumsAtom1[Wrapper.t1Rhoj0Bin].PE(t1GradRhoj0);
    	
    	//t2GradRhoj0
    	t2GradRhoj0.Ea1Tv1(p.t2, gradRhoj0);
    	
    	gradientSumsAtom0[Wrapper.t2Rhoj0Bin].PE(t2GradRhoj0);
    	gradientSumsAtom1[Wrapper.t2Rhoj0Bin].PE(t2GradRhoj0);
    	
    	//t3GradRhoj0
    	t3GradRhoj0.Ea1Tv1(p.t3, gradRhoj0);
    	
    	gradientSumsAtom0[Wrapper.t3Rhoj0Bin].PE(t3GradRhoj0);
    	gradientSumsAtom1[Wrapper.t3Rhoj0Bin].PE(t3GradRhoj0);
    	
    	
    	
    	//Calculations to determine gradient of phi
    	
    	gradRhoi0Ref.Ea1Tv1(p.Z, gradRhoj0);
    	
    	gradRhoi1Ref.Ea1Tv1(x, gradRhoj1x);
    	gradRhoi1Ref.PEa1Tv1(y, gradRhoj1y);
    	gradRhoi1Ref.PEa1Tv1(z, gradRhoj1z);
    	gradRhoi1Ref.TE(p.Z*p.Z*rhoj1/rhoi1Ref);
    	
    	gradRhoi2Ref.Ea1Tv1(x*x, gradRhoj2xx);
    	gradRhoi2Ref.PEa1Tv1(2*x*y, gradRhoj2xy);
    	gradRhoi2Ref.PEa1Tv1(2*x*z, gradRhoj2xz);
    	gradRhoi2Ref.PEa1Tv1(y*y, gradRhoj2yy);
    	gradRhoi2Ref.PEa1Tv1(2*y*z, gradRhoj2yz);
    	gradRhoi2Ref.PEa1Tv1(z*z, gradRhoj2zz);
    	gradRhoj213.Ea1Tv1(1/3, gradRhoj2);
    	gradRhoi2Ref.ME(gradRhoj213);
    	gradRhoi2Ref.TE(p.Z*p.Z*rhoj2/rhoi2Ref);
    	
    	gradRhoi3Ref.Ea1Tv1(x*x*x, gradRhoj3xxx);
    	gradRhoi3Ref.PEa1Tv1(3*x*x*y, gradRhoj3xxy);
    	gradRhoi3Ref.PEa1Tv1(3*x*x*z, gradRhoj3xxz);
    	gradRhoi3Ref.PEa1Tv1(3*x*y*y, gradRhoj3xyy);
    	gradRhoi3Ref.PEa1Tv1(6*x*y*z, gradRhoj3xyz);
    	gradRhoi3Ref.PEa1Tv1(3*x*z*z, gradRhoj3xzz);
    	gradRhoi3Ref.PEa1Tv1(y*y*y, gradRhoj3yyy);
    	gradRhoi3Ref.PEa1Tv1(3*y*y*z, gradRhoj3yyz);
    	gradRhoi3Ref.PEa1Tv1(3*y*z*z, gradRhoj3yzz);
    	gradRhoi3Ref.PEa1Tv1(z*z*z, gradRhoj3zzz);
    	gradRhoi3Ref.TE(p.Z*p.Z*rhoj3/rhoi3Ref);
    	
    	gradGammaRef.Ea1Tv1( (-1/rhoi0Ref) * (p.t1*rhoi1Ref*rhoi1Ref +
    			p.t2*rhoi2Ref*rhoi2Ref + p.t3*rhoi3Ref*rhoi3Ref), gradRhoi0Ref);
    	gradGammaRef.PEa1Tv1(p.t1*rhoi1Ref, gradRhoi1Ref);
    	gradGammaRef.PEa1Tv1(p.t2*rhoi2Ref, gradRhoi2Ref);
    	gradGammaRef.PEa1Tv1(p.t3*rhoi3Ref, gradRhoi3Ref);
    	gradGammaRef.TE(2/(rhoi0Ref*rhoi0Ref));
    	
    	gradRhoiRef.Ea1Tv1(rhoi0Ref*Math.exp(-gammaRef)
    			/((1+Math.exp(-gammaRef))*(1+Math.exp(-gammaRef))), gradGammaRef);
    	gradRhoiRef.Ea1Tv1(1/(1+Math.exp(-gammaRef)), gradRhoi0Ref);
    	gradRhoiRef.TE(2);
    	
    	gradFRef.Ea1Tv1( (p.A*p.Ec/p.rhoScale)*
    			(1+Math.log(rhoiRef)-Math.log(p.rhoScale)), gradRhoiRef);
    	
    	gradERef.Ea1Tv1( (p.Ec*p.alpha*p.alpha/p.r0) * ((r/p.r0) - 1)
    					*(Math.exp(-p.alpha*((r/p.r0)-1))), gradr);
    	
    	gradPhi.E(gradERef);
    	gradPhi.ME(gradFRef);
    	gradPhi.TE(2/p.Z);
    	
    	gradientSumsAtom0[Wrapper.phiBin].PE(gradPhi);
    	gradientSumsAtom1[Wrapper.phiBin].PE(gradPhi);
    	
    	//MEAMP2 should contribute nothing to the gradient directly.
        return nada;
    }
    
    /**
     * Same as uInt.
     */
    public double integral(double rC) {
        return 0;
    }
    
    /**
     * Returns infinity.  May be overridden to define a finite-ranged potential.
     */
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public PhaseAgentManager getPhaseAgentManager(){
    	return phaseAgentManager;
    }
    
    private final Vector3D unitVector;
    private final Vector3D vector100;
    private final Vector3D vector010;
    private final Vector3D vector001;
    private final Vector3D gradx;
    private final Vector3D grady;
    private final Vector3D gradz;
    private final Vector3D gradr;
    private final Vector3D gradRhoj0;
    private final Vector3D gradRhoj1;
    private final Vector3D gradRhoj2;
    private final Vector3D gradRhoj3;
    private final Vector3D gradRhoj1x;
    private final Vector3D gradRhoj1y;
    private final Vector3D gradRhoj1z;
    private final Vector3D gradRhoj2xx;
    private final Vector3D gradRhoj2xy;
    private final Vector3D gradRhoj2xz;
    private final Vector3D gradRhoj2yy;
    private final Vector3D gradRhoj2yz;
    private final Vector3D gradRhoj2zz;
    private final Vector3D gradRhoj3xxx;
    private final Vector3D gradRhoj3xxy;
    private final Vector3D gradRhoj3xxz;
    private final Vector3D gradRhoj3xyy;
    private final Vector3D gradRhoj3xyz;
    private final Vector3D gradRhoj3xzz;
    private final Vector3D gradRhoj3yyy;
    private final Vector3D gradRhoj3yyz;
    private final Vector3D gradRhoj3yzz;
    private final Vector3D gradRhoj3zzz;
    private final Vector3D t1GradRhoj0;
    private final Vector3D t2GradRhoj0;
    private final Vector3D t3GradRhoj0;
    private final Vector3D gradRhoi0Ref;
    private final Vector3D gradRhoi1Ref;
    private final Vector3D gradRhoi2Ref;
    private final Vector3D gradRhoj213;
    private final Vector3D gradRhoi3Ref;
    private final Vector3D gradGammaRef;
    private final Vector3D gradRhoiRef;
    private final Vector3D gradERef;
    private final Vector3D gradFRef;
    private final Vector3D gradPhi;
    private final Vector3D nada;
    
    private final CoordinatePair coordinatePair;
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
    private ParameterSetMEAM p;
    private double r2Last = -1.0;
    private AtomPair pair;
    protected AtomAgentManager[] agentManager;
    protected Wrapper[] agents;
    private final PhaseAgentManager phaseAgentManager;
    
    public Class getAgentClass() {
        return Wrapper.class;
    }
    
    public Object makeAgent(Atom atom) {
        return new Wrapper();
    }
    
    public void releaseAgent(Object agent, Atom atom) {}
    
    public static class Wrapper implements java.io.Serializable {
    	public Wrapper() {
    		sums = new double[27];
    		gradientSums = new Vector3D[27];
    		for (int i=0; i<27; i++){
    			gradientSums[i] = new Vector3D();
    		}
    	}
    	public double[] sums;
    	public Vector3D[] gradientSums;
    	public static int rhoj0Bin = 0;
    	public static int t1Rhoj0Bin = 1;
    	public static int t2Rhoj0Bin = 2;
    	public static int t3Rhoj0Bin = 3;
    	public static int rhoj1Bin = 4;
    	public static int rhoj1xBin = 5;
    	public static int rhoj1yBin = 6;
    	public static int rhoj1zBin = 7;
    	public static int rhoj2Bin = 8;
    	public static int rhoj2xxBin = 9;
    	public static int rhoj2xyBin = 10;
    	public static int rhoj2xzBin = 11;
    	public static int rhoj2yyBin = 12;
    	public static int rhoj2yzBin = 13;
    	public static int rhoj2zzBin = 14;
    	public static int rhoj3Bin = 15;
    	public static int rhoj3xxxBin = 16;
    	public static int rhoj3xxyBin = 17;
    	public static int rhoj3xxzBin = 18;
    	public static int rhoj3xyyBin = 19;
    	public static int rhoj3xyzBin = 20;
    	public static int rhoj3xzzBin = 21;
    	public static int rhoj3yyyBin = 22;
    	public static int rhoj3yyzBin = 23;
    	public static int rhoj3yzzBin = 24;
    	public static int rhoj3zzzBin = 25;
    	public static int phiBin = 26;
    }
    
}
