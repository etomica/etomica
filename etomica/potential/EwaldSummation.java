package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomLeafAgentManager;
import etomica.math.SpecialFunctions;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.units.Electron;
//import etomica.potential.EwaldSumMolecules2.MyCharge;
/**
 * basic Ewald Sum 
 * U(coulomb) = U(real-space) + U(fourier-space) + U(self-correction)
 * Frenkel&Smit, <Understanding Molecular Simulation from Algorithms to Applications> 
 * p292-300.especially equation(12.1.25)
 * also Allen, computer simulation of liquids Ch5
 * 
 * 
 * @author shu
 * Date:08-09-2012
 */

/**Sabry*/
//Given: rc=L/2  ,  s  ,  
//Calc.: alpha = s/rc    ,  F&S: radius => nc=(s*alpha/PI)*L Sabry: n_max=N  , nx = ny = nz ~ N^1/3   ,
//                            Shu : nx = ny = nz = coefficient_fourier * N^1/3 (coefficient_fourier=4 hERE)

//OLD: public class EwaldSumMolecules implements IPotentialMolecular {
public class EwaldSummation implements PotentialSoft{
		
	public final ISpace space;
	public final AtomLeafAgentManager atomAgentManager;
	public final IBox box;
	public final double alpha;//sqrt of the Frenkel's alpha, follow other authors' convention
	public final double boxSize;
	public final double volume;
	public final double rCut; // real-space spherical cutoff
	
	public final double precision_s;
	public final double exp_s;// = exp(-s*s)/s/s;
	public final int coefficient_fourier = 4; 
	public final IMoleculeList moleculeList;
	public final int numMolecules;
	public final double q_err; // numberMolecules * charge in sim unit, for err estimate
	protected IVectorMutable[] gradient;
	protected double[] sinkrj, coskrj;
	protected int nc;
		
	// *********************************************** constructor ************************************ // 
	public EwaldSummation(IBox box, AtomLeafAgentManager atomAgentManager, double precision_s, ISpace _space){
		this.box = box;
		this.atomAgentManager = atomAgentManager;
		this.precision_s = precision_s; // =???? | SG
		this.space = _space;
		boxSize = box.getBoundary().getBoxSize().getX(0);
		volume = box.getBoundary().volume();
		double foo = 2.0;
		rCut = boxSize / foo ; 
		alpha = precision_s / rCut;//Separation parameter, obtained from s and L
		nc = 20; // Converged enough for sI :) 
				
		double precision_sSquared = precision_s * precision_s ;
		exp_s = Math.exp(-precision_sSquared) / precision_sSquared;
		moleculeList = box.getMoleculeList();
		numMolecules = moleculeList.getMoleculeCount();
		//nc = (int)Math.ceil(precision_s * boxSize * alpha / Math.PI);//"n = 2s^2/pi" if rC = L/2

		double e = Electron.UNIT.toSim(1.0); 
		q_err = numMolecules * e * e ;
		
//		System.out.println("In Ewald sum class");
//		System.out.println("exp(-s^2)/s/s:"+exp_s);
//		System.out.println("box size is : "+boxSize);
//		System.out.println("rCut = L /"+foo+":"+rCut);
//		System.out.println("alpha = s / rCut:"+alpha);
//		System.out.println("cf. 5/L:"+5.0/boxSize);
//		System.out.println("    3.5/L:"+3.5/boxSize);
//		System.out.println("q_err= N*q^2:"+ q_err);
		gradient = new IVectorMutable[0];
		sinkrj = new double[0];
		coskrj = new double[0];
	}
	
	//////////////////////////////////////////// begin calculating energy //////////////////////////////////////
	
	// *********************************************************************************************//
	// *************************************  Real-space ******************************************//
	// *********************************************************************************************//
	public double uReal(){
		int nAtoms = box.getLeafList().getAtomCount();
		
		double rCutSquared = rCut * rCut; // criteria for spherical cutoff
		double uReal = 0.0;
		
		IVectorMutable rAB = space.makeVector();// vector between site A @ molecule i & site B @ molecule j
		for (int i=0; i < nAtoms; i++){//H
			IAtom atomA = box.getLeafList().getAtom(i);
			int aIndex = atomA.getParentGroup().getIndex();
			IVectorMutable positionA = atomA.getPosition();
			double chargeA = ((MyCharge)atomAgentManager.getAgent(atomA)).charge;
			for (int j=i+1; j < nAtoms; j++){
				IAtom atomB = box.getLeafList().getAtom(j);
				int bIndex = atomB.getParentGroup().getIndex();
				
				if(aIndex == bIndex) continue;//Skip same molecules!
				
				IVectorMutable positionB = atomB.getPosition();
				double chargeB = ((MyCharge)atomAgentManager.getAgent(atomB)).charge;
				
 				rAB.Ev1Mv2(positionA, positionB);// get vector rAB
				box.getBoundary().nearestImage(rAB);// minimum image
				double rABMagnitudeSquared = rAB.squared() ;// Squared of |rAB|
				if (rABMagnitudeSquared < 0.5){
					throw new RuntimeException();
				}
				if (rABMagnitudeSquared > rCutSquared) continue ; // check whether the rAB is within the spherical cutoff
				double rABMagnitude = Math.sqrt(rABMagnitudeSquared);

				uReal += chargeA * chargeB * SpecialFunctions.erfc( alpha * rABMagnitude) / rABMagnitude;//Don't worry about 1/2 factor!
			}// close for all sites in j-th molecule
		} // close for the outside loop
		return uReal;
	}

	// *********************************************************************************************//
	// *************************************  Fourier-space ****************************************//
	// *********************************************************************************************//
	public double uFourier(){
//		System.out.println("coefficient_fourier: " + coefficient_fourier); 
//		System.out.println("number of k-vectors along one axis is :"+ nc); 
		// int n =(int)( precision_s * boxSize * alpha / Math.PI);// n(cut): number of vectors along one axis in fourier-space
		
		double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
		double kCut =  basis * nc;
		double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
		double coefficient = 2.0 / volume * Math.PI ;
		//double uFourier = 0.0;
		double uFourier = 0.0; //>>>>>>>>>>>>>> calculated from cos*cos + sin*sin
		IVectorMutable kVector = space.makeVector();// fourier space vector
		IVectorMutable r = space.makeVector();
		int nAtoms = box.getLeafList().getAtomCount();

		if (!true){
			System.out.println("basis vector is :" + basis); 
			System.out.println("kCutoff is :" + kCut); 
		}
		
		// loop over vectors in k-space. (1)k=(0,0,0) is excluded, (2)within sphere with kCutoff as its radius
		for (int xAxis = -nc; xAxis < nc+1; xAxis++){
			kVector.setX(0, (xAxis * basis));// assign value to the x-axis
			for (int yAxis = -nc; yAxis < nc+1; yAxis++ ){
				kVector.setX(1, (yAxis * basis));// assign value to the y-axis
				for (int zAxis = -nc; zAxis < nc+1; zAxis++ ){
					if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
					kVector.setX(2, (zAxis * basis));// assign value to the z-axis, now the vector is specified
					double kSquared = kVector.squared();
					
					if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
					
					double kCoefficientTerm = Math.exp(-kSquared / 4.0 / alpha / alpha) / kSquared;//exp(-k*k/4/alpha/alpha) / k/k , a constant for a given kVector
					double structureFactorReal =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
					double structureFactorImagine =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
					for (int i=0; i< nAtoms ; i++){
						IAtom atom = box.getLeafList().getAtom(i);
						IVectorMutable position = atom.getPosition();
						double charge = ((MyCharge)atomAgentManager.getAgent(atom)).charge;
							r.E(position);
							double k_r_dot = kVector.dot(r);
							structureFactorReal += charge * Math.cos(k_r_dot);// >>>>>>>>>>>>> calculated from cos*cos + sin*sin
							structureFactorImagine += charge * Math.sin(k_r_dot);// >>>>>>>>>>>>> calculated from cos*cos + sin*sin
					}// close loop for all sites within 1 molecule
					double structureFactorSquared = structureFactorReal * structureFactorReal + structureFactorImagine * structureFactorImagine ;//>>>>>>>>>>>>>calculated from cos*cos + sin*sin
					uFourier +=  kCoefficientTerm * structureFactorSquared; //>>>>>>>>>>>>>calculated from cos*cos + sin*sin
				}// close for z-axis
			}//close for y-axis
		}// close for x-axis(all non-zero kVecors)
		double u = coefficient * uFourier; 

		return u; 
	}
	
	
	// *********************************************************************************************//
	// ********************** self-correction Part************************************************* // 
	// *********************************************************************************************//
	public double uSelf(){
		double coefficient = -alpha/Math.sqrt(Math.PI);
		double uSelf = 0.0;

		boolean RPM1_1_78point5 = false;
		if (RPM1_1_78point5 ){
//			System.out.println("dielectric permitivity is 78.5,and the system is 1-1 RPM, so I am using u_self = -alpha/sqrt(pi) * number * e*e");
			double eCharge =  Electron.UNIT.toSim( 1.0/Math.sqrt(78.5)); 
			uSelf = coefficient * numMolecules * eCharge * eCharge ; 
		} 
		else {
//			System.out.println("More generic algorithm, and I am looping every ion in the system");
			for (int i=0; i< numMolecules; i++){
				IMolecule molecule = moleculeList.getMolecule(i);	
				int numSites = molecule.getChildList().getAtomCount();
				// each site has a charge. get charge info from every site, loop over all sites
				for (int site=0; site<numSites; site++){
					IAtom atom = molecule.getChildList().getAtom(site);					
					double charge = ((MyCharge)atomAgentManager.getAgent(atom)).charge;
//					uSelf += coefficient*charge*charge;
					uSelf += charge*charge;

				}
			}
		}
		
		uSelf *= coefficient;
		return uSelf;
	}
	public double uBondCorr(){
		double uCorr = 0.0;
		IVectorMutable rAB = space.makeVector();// vector between site A @ molecule i & site B @ molecule j
		for (int i=0; i< numMolecules; i++){
			IMolecule molecule = moleculeList.getMolecule(i);	
			int numSites = molecule.getChildList().getAtomCount();
			for (int siteA=0; siteA<numSites; siteA++){
				IAtom atomA = molecule.getChildList().getAtom(siteA);
				IVectorMutable positionA = atomA.getPosition();
				double chargeA = ((MyCharge)atomAgentManager.getAgent(atomA)).charge;
				for (int siteB=siteA+1; siteB<numSites; siteB++){
					IAtom atomB = molecule.getChildList().getAtom(siteB);
					IVectorMutable positionB = atomB.getPosition();
					double chargeB = ((MyCharge)atomAgentManager.getAgent(atomB)).charge;

	 				rAB.Ev1Mv2(positionA, positionB);
					box.getBoundary().nearestImage(rAB);
					double rABMagnitudeSquared = rAB.squared();
					double rABMagnitude = Math.sqrt(rABMagnitudeSquared);
					
					uCorr += chargeA*chargeB*(1-SpecialFunctions.erfc(alpha*rABMagnitude))/rABMagnitude;
				}
			}
		}		
		return uCorr;
	}
	// ************************************   non-reduced energy   ******************************************************
	public double sum(){
		if (false){
			long t1_real = System.currentTimeMillis();
			for ( int i = 0 ;i < 100000; i++){
				double real = uReal();
			}
			long t2_real = System.currentTimeMillis();
//			System.out.println("delta_t_real: "+ (t2_real- t1_real));
		}
		if (false){
			long t1_fourier = System.currentTimeMillis();
			for ( int i = 0 ;i < 10000; i++){
				double fourier = uFourier();
			}
			long t2_fourier = System.currentTimeMillis();
//			System.out.println("delta_t_fourier: "+ (t2_fourier- t1_fourier));
		}
		
		double real = uReal();
		double fourier = uFourier();
		double self = uSelf();
		double bondCorr = uBondCorr();

		double nonReducedEnergy = real + fourier + self - bondCorr;
			
		if (!false) { 
//			System.out.println("total:               "+ nonReducedEnergy/numMolecules);
//			System.out.println("real   : "+ real/numMolecules);
//			System.out.println("fourier: "+ fourier/numMolecules);
//			System.out.println("self: "+ self/numMolecules);
		}
		
		return nonReducedEnergy;
	}
	
	
	public double energy(IMoleculeList atoms) {
		double energy = sum();	
		//System.out.print("Total:  "+ energy);
		return energy;
	}


	
	public double getRange() {
		// TODO Auto-generated method stub
		return Double.POSITIVE_INFINITY;
	}

	public int nBody() {
		// TODO Auto-generated method stub
		return 0;
	}

	public void setBox(IBox box) {
		// do nothing
		
	}
	//******************************** inner class ********************************************//
	public static class MyCharge{
		public MyCharge(double charge){
			this.charge = charge;
		}
		public final double charge;
		
	}
	
	@Override
	public double energy(IAtomList atoms) {
		// TODO Auto-generated method stub
		double energy = sum();	
		//System.out.print("Total:  "+ energy);
		return energy;
	}

	@Override
	public double virial(IAtomList atoms) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
	//////////////////////////////////////////// begin calculating gradient //////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	
	@Override
	public IVector[] gradient(IAtomList atoms) {
		// TODO Auto-generated method stub
		int nAtoms = box.getLeafList().getAtomCount();
		double rCutSquared = rCut * rCut; // criteria for spherical cutoff
		IVectorMutable rAB = space.makeVector();// vector between site A @ molecule i & site B @ molecule j
		double coeff = 4.0*Math.PI / box.getBoundary().volume();
		double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
		double kCut =  basis * nc;
		double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
		IVectorMutable kVector = space.makeVector();// fourier space vector

		if(gradient.length < nAtoms){
			gradient = new IVectorMutable[nAtoms];
			sinkrj = new double[nAtoms];
			coskrj = new double[nAtoms];
			for(int i=0; i < nAtoms; i++){
				gradient[i] = space.makeVector();	
			}
		}else{
			for(int i=0; i < nAtoms; i++){
				gradient[i].E(0); //Vecors and scalar	
				sinkrj[i] = 0;
				coskrj[i] = 0;
			}	
		}
		
//Real gradient  //Cross Interaction
		for (int i=0; i < nAtoms; i++){
			IAtom atomA = box.getLeafList().getAtom(i);
			int aIndex = atomA.getParentGroup().getIndex(); // molecule a
			IVectorMutable positionA = atomA.getPosition();
			double chargeA = ((MyCharge)atomAgentManager.getAgent(atomA)).charge;
				for (int j=i+1; j < nAtoms; j++){
					IAtom atomB = box.getLeafList().getAtom(j);
					int bIndex = atomB.getParentGroup().getIndex(); // molecule b
					
					if(aIndex == bIndex) continue;//Skip same molecules!

					IVectorMutable positionB = atomB.getPosition();
					rAB.Ev1Mv2(positionA, positionB); //rAB == rA - rB
					box.getBoundary().nearestImage(rAB);
					double rAB2 = rAB.squared() ;
					if (rAB2 > rCutSquared) continue ; 
					double chargeB = ((MyCharge)atomAgentManager.getAgent(atomB)).charge;
					double rABMagnitude = Math.sqrt(rAB2);
					double rAB3 = rABMagnitude*rAB2;
					double B = SpecialFunctions.erfc(alpha*rABMagnitude) + 2.0*alpha*rABMagnitude/Math.sqrt(Math.PI) * Math.exp(-alpha*alpha*rAB2) ;
					double realCoeff = - chargeA*chargeB * B / rAB3; // gradU = -F
			        gradient[i].PEa1Tv1(realCoeff, rAB);
			        gradient[j].PEa1Tv1(-realCoeff, rAB);
//			        if(i==0){
//			            System.out.println("Real = " + gradient[0]);	
//			        }
			        
				}
			}

		
		
//Fourier gradient Part		

		for (int xAxis = -nc; xAxis < nc+1; xAxis++){
			kVector.setX(0, (xAxis * basis));// assign value to the x-axis
			for (int yAxis = -nc; yAxis < nc+1; yAxis++ ){
				kVector.setX(1, (yAxis * basis));// assign value to the y-axis
				for (int zAxis = -nc; zAxis < nc+1; zAxis++ ){
					if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
					kVector.setX(2, (zAxis * basis));// assign value to the z-axis, now the vector is specified
					double kSquared = kVector.squared();
					if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
					double sCoskr = 0.0, sSinkr = 0.0;
					
					for (int j=0; j< nAtoms ; j++){ // Loop over atoms (4*nBasis)
						IAtom atom = box.getLeafList().getAtom(j);
						IVectorMutable position = atom.getPosition();
						sinkrj[j] = Math.sin(position.dot(kVector));
						coskrj[j] = Math.cos(position.dot(kVector));
						double chargej = ((MyCharge)atomAgentManager.getAgent(atom)).charge;

						sCoskr += chargej*coskrj[j];
						sSinkr += chargej*sinkrj[j]; 
					}//End loop over j
					
					double coeffk = coeff / kSquared * Math.exp(-kSquared/4.0/alpha/alpha);
					for(int i=0; i<nAtoms; i++){
						IAtom atom = box.getLeafList().getAtom(i);
                		double chargei = ((MyCharge)atomAgentManager.getAgent(atom)).charge;
						double coeffki = coeffk * chargei * (sinkrj[i] * sCoskr - coskrj[i] * sSinkr); 
						gradient[i].PEa1Tv1(-coeffki , kVector);  // gradU = -F
					}
//		            System.out.println("Fourier = " + gradient[0]);	

				}//end of storing Sin and Cos
			}
		}//End loop over ks
		
		
//Intra-Molecular  gradient:
		for (int i=0; i< numMolecules; i++){
			IMolecule molecule = moleculeList.getMolecule(i);	
			int numSites = molecule.getChildList().getAtomCount();
			for (int siteA=0; siteA<numSites; siteA++){
				IAtom atomA = molecule.getChildList().getAtom(siteA); // index = 0, 1, 2, 3|||leafIndex=0...184
				IVectorMutable positionA = atomA.getPosition();
				double chargeA = ((MyCharge)atomAgentManager.getAgent(atomA)).charge;
				for (int siteB=siteA+1; siteB<numSites; siteB++){
					IAtom atomB = molecule.getChildList().getAtom(siteB);
					IVectorMutable positionB = atomB.getPosition();
					double chargeB = ((MyCharge)atomAgentManager.getAgent(atomB)).charge;

	 				rAB.Ev1Mv2(positionA, positionB);
					box.getBoundary().nearestImage(rAB);
					double rAB2 = rAB.squared();
					double rABMagnitude = Math.sqrt(rAB2);
					// U = Ur + Uf - Uself - U_intra ====> U_intra = Erf(alpha r)/r
					// dU = -F = Ur' + Uf' - Uself' - U_intra' ===>   -d[Erf(alpha r)/r]/dx
					double B = 2*alpha/Math.sqrt(Math.PI) * Math.exp(-alpha*alpha*rAB2)-(1-SpecialFunctions.erfc(alpha*rABMagnitude))/rABMagnitude; 
					double coeffAB = - chargeA*chargeB * B / rAB2; // gradU = -F
//					System.out.println(atomA.getLeafIndex());
//					System.out.println(atomB.getLeafIndex());
			        gradient[atomA.getLeafIndex()].PEa1Tv1(coeffAB, rAB);
			        gradient[atomB.getLeafIndex()].PEa1Tv1(-coeffAB, rAB);
//			        if(i==0 && siteA == 0){
//			            System.out.println("Intra = " + gradient[0]);	
//			        }

				}
			}
		}		
		return gradient;
	}

	@Override
	public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
		// TODO Auto-generated method stub
		return null;
	}
	
	
}
