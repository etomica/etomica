/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.*;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager;
import etomica.math.Complex;
import etomica.math.SpecialFunctions;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.Electron;

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

public class EwaldSumMolecules implements IPotentialMolecular {
		
	public final Space space;
	public final AtomLeafAgentManager atomAgentManager;
	public final Box box;
	public final double temperature;// in simulation unit
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
		
	// *********************************************** constructor ************************************ // 
	public EwaldSumMolecules(Box box, AtomLeafAgentManager atomAgentManager, double precision_s, double temperature, Space _space){
		this.box = box;
		this.atomAgentManager = atomAgentManager;
		this.precision_s = precision_s;
		this.space = _space;
		this.temperature = temperature;
		
		boxSize = box.getBoundary().getBoxSize().getX(0);
		volume = box.getBoundary().volume();
		double foo = 2.0;
		rCut = boxSize / foo ; 
		alpha = precision_s / rCut;//Separation parameter, obtained from s and L
		double precision_sSquared = precision_s * precision_s ;
		exp_s = Math.exp(-precision_sSquared) / precision_sSquared;
		moleculeList = box.getMoleculeList();
		numMolecules = moleculeList.getMoleculeCount();

		double e = Electron.UNIT.toSim(1.0); 
		q_err = numMolecules * e * e ;
		
		System.out.println("In Ewald sum class");
		System.out.println("exp(-s^2)/s/s:"+exp_s);
		System.out.println("box size is : "+boxSize);
		System.out.println("rCut = L /"+foo+":"+rCut);
		System.out.println("alpha = s / rCut:"+alpha);
		System.out.println("cf. 5/L:"+5.0/boxSize);
		System.out.println("    3.5/L:"+3.5/boxSize);
		System.out.println("q_err= N*q^2:"+ q_err);
		
	}
	
	//////////////////////////////////////////// begin calculating energy //////////////////////////////////////
	
	// *********************************************************************************************//
	// *************************************  Real-space ******************************************//
	// *********************************************************************************************//
	public double uReal(){
		
		double rCutSquared = rCut * rCut; // criteria for spherical cutoff
		double uReal = 0.0;
		
		Vector rAB = space.makeVector();// vector between site A @ molecule i & site B @ molecule j
		for (int i=0; i < numMolecules; i++){
			IMolecule molecule_i = moleculeList.getMolecule(i); // get i-th molecule
			int numSites = molecule_i.getChildList().getAtomCount();
			
			for (int a=0; a < numSites; a++){
				IAtom siteA = molecule_i.getChildList().getAtom(a);// get siteA from i-th molecule
				Vector positionA = siteA.getPosition();
				double chargeA = ((MyCharge)atomAgentManager.getAgent(siteA)).charge;
				
				// given i-th molecule, get j-th molecule starting from (i+1)-th molecule
				for (int j=i+1; j < numMolecules; j++){
					
					IMolecule molecule_j = moleculeList.getMolecule(j);
					// get siteB from molecule_j
					for (int b=0; b < numSites ; b++){
						
						IAtom siteB = molecule_j.getChildList().getAtom(b);
						Vector positionB = siteB.getPosition();
						double chargeB = ((MyCharge)atomAgentManager.getAgent(siteB)).charge;
						
						rAB.Ev1Mv2(positionA, positionB);// get vector rAB
						box.getBoundary().nearestImage(rAB);// minimum image
						
						double rABMagnitudeSquared = rAB.squared() ;// Squared of |rAB|
						if (rABMagnitudeSquared > rCutSquared) continue ; // check whether the rAB is within the spherical cutoff
						double rABMagnitude = Math.sqrt(rABMagnitudeSquared);
						uReal += chargeA * chargeB * SpecialFunctions.erfc( alpha * rABMagnitude) / rABMagnitude;//Don't worry about 1/2 factor!
					}// close for all sites in j-th molecule
				
				}//close for all molecule j for a given siteA in i-th molecule
			}// close for all sites in i-th molecule
		} // close for the outside loop
		
		// ************************************** error in real-space calculation ***************************************************
		// *****************From Kolafa & Perram, cutoff errors in the Ewald summation formulae for point charge systems************
		//double errReal =  q_err * Math.sqrt(precision_s / 2 / alpha / volume) * exp_s; 
		//System.out.println("*********************Real-space ***********"); 
		//System.out.println("error_uReal:"+errReal); 
		//System.out.println("relative error of real-space energy calculation:"+errReal / uReal); 
//		System.out.println("uReal:             "+uReal); 

		return uReal;
	}
	

	// *********************************************************************************************//
	// *************************************  fourier-space ****************************************//
	// *********************************************************************************************//
	public double uFourier(){
		int n0 = (int)Math.round(Math.pow(numMolecules / 2.0, 1.0/3) ) ;
		int n = n0 * coefficient_fourier; 
		System.out.println("number of vectors along one real space axis is: "+n0); 
		System.out.println("coefficient_fourier: " + coefficient_fourier); 
		System.out.println("number of k-vectors along one axis is :"+ n); 
		// int n =(int)( precision_s * boxSize * alpha / Math.PI);// n(cut): number of vectors along one axis in fourier-space
		
		double basis = 2 * Math.PI / boxSize; // basis(unit) vector magnitude of fourier space vector
		double kCut =  basis * n;
		double kCutSquared = kCut * kCut; // criteria for spherical cutoff in fourier space
		double coefficient = 2.0 / volume * Math.PI ;
		boolean useComplex = true;
		double uFourier = 0.0;
		double uFourier_ = 0.0; //>>>>>>>>>>>>>> calculated from cos*cos + sin*sin
		Vector kVector = space.makeVector();// fourier space vector
		Vector r = space.makeVector();
		
		if (true){
			System.out.println("basis vector is :" + basis); 
			System.out.println("kCutoff is :" + kCut); 
		}
		
		// loop over vectors in k-space. (1)k=(0,0,0) is excluded, (2)within sphere with kCutoff as its radius
		for (int xAxis = -n; xAxis < n+1; xAxis++){
			kVector.setX(0, (xAxis * basis));// assign value to the x-axis
			
			for (int yAxis = -n; yAxis < n+1; yAxis++ ){
				kVector.setX(1, (yAxis * basis));// assign value to the y-axis
				
				for (int zAxis = -n; zAxis < n+1; zAxis++ ){
					
					if( (xAxis * xAxis + yAxis * yAxis + zAxis * zAxis) == 0) continue;// first check: k is a non-zero vector
					
					kVector.setX(2, (zAxis * basis));// assign value to the z-axis, now the vector is specified
					
					double kSquared = kVector.squared();
					if (kSquared > kCutSquared) continue;// k-vector should be within the sphere with kCutoff as the radius
					
					double kCoefficientTerm = Math.exp(-kSquared / 4.0 / alpha / alpha) / kSquared;//exp(-k*k/4/alpha/alpha) / k/k , a constant for a given kVector
					
					/* calculate |rho(k)|^2, structureFactor
					 * for a given kVector, loop over all rVectors involved in real-space
					 * structureFactor = sum [ qi * exp(-i * kVectors <dot product> rVectors) ]     
					 */
					Complex structureFactor = new Complex(0.0, 0.0);
					double structureFactorReal =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
					double structureFactorImagine =  0.0;//>>>>>>>>>>>>> calculated from cos*cos + sin*sin
					
					// get molecule_i
					for (int i=0; i< numMolecules ; i++){
						
						Complex expInside = new Complex(0.0, 0.0); 
						Complex expWithCharge = new Complex(0.0, 0.0);
						
						IMolecule molecule_i = moleculeList.getMolecule(i);
						int numSites = molecule_i.getChildList().getAtomCount();
						
						// get the interaction site from molecule_i
						for (int a=0; a< numSites ; a++){
							IAtom site = molecule_i.getChildList().getAtom(a);
							Vector position = site.getPosition();
							double charge = ((MyCharge)atomAgentManager.getAgent(site)).charge;
							r.E(position);
							double k_r_dot = kVector.dot(r);
							
							structureFactorReal += charge * Math.cos(k_r_dot);// >>>>>>>>>>>>> calculated from cos*cos + sin*sin
							structureFactorImagine += charge * Math.sin(k_r_dot);// >>>>>>>>>>>>> calculated from cos*cos + sin*sin
							
							expInside = new Complex(0.0,k_r_dot);
							expWithCharge = expInside.exponential().times(new Complex(charge,0.0));
							structureFactor = structureFactor.plus(expWithCharge); // i.e. strutureFactor += expWithCharge
							
						}// close loop for all sites within 1 molecule
						
					}// close for all molecules corresponding to 1 kVector
					double modStructureFactor = structureFactor.modulus();// for a given k-vector, sum of all r-vectors
					double modStructureFactorSquared = modStructureFactor * modStructureFactor; 
					uFourier += ( kCoefficientTerm * modStructureFactorSquared ) ; 
					
					double structureFactorSquared = structureFactorReal * structureFactorReal + structureFactorImagine * structureFactorImagine ;//>>>>>>>>>>>>>calculated from cos*cos + sin*sin
					uFourier_ +=  kCoefficientTerm * structureFactorSquared; //>>>>>>>>>>>>>calculated from cos*cos + sin*sin
				}// close for z-axis
			}//close for y-axis
		}// close for x-axis(all non-zero kVecors)
		double u = coefficient * uFourier; 
		// ************************************** errors in fourier space ***********************
		//double errFourier = q_err  *  Math.sqrt(precision_s / Math.PI/ alpha / volume)  * exp_s; 
//		System.out.println("uFourier:           "+u); 
		//System.out.println("error_fourier:"+errFourier); 
		//System.out.println("relative error of fourier-space energy calculation:"+errFourier / u); 
		return u; 
	}
	
	// *********************************************************************************************//
	// ********************** self-correction Part************************************************* // 
	// *********************************************************************************************//
	// U(self)= [-alpha/sqrt(pi)]* sum of (qi*qi), i: all ions in the system
	// U(self) is a constant for a given system regardless of the configurations
	public double uSelf(){
		double coefficient = -alpha / Math.sqrt(Math.PI);
		double uSelf = 0.0;

		boolean RPM1_1_78point5 = false;
		if (RPM1_1_78point5 ){
			System.out.println("dielectric permitivity is 78.5,and the system is 1-1 RPM, so I am using u_self = -alpha/sqrt(pi) * number * e*e");
			double eCharge =  Electron.UNIT.toSim( 1.0/Math.sqrt(78.5)); 
			uSelf = coefficient * numMolecules * eCharge * eCharge ; 
		} 
		else {
			System.out.println("More generic algorithm, and I am looping every ion in the system");
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
		
//		System.out.println("uSelf:             " + uSelf);
//		return uSelf;
		uSelf *= coefficient;
		return uSelf;

	}
	// ************************************   non-reduced energy   ******************************************************
	public double sum(){
		if (false){
			long t1_real = System.currentTimeMillis();
			for ( int i = 0 ;i < 100000; i++){
				double real = uReal();
			}
			long t2_real = System.currentTimeMillis();
			System.out.println("delta_t_real: "+ (t2_real- t1_real));
		}
		if (false){
			long t1_fourier = System.currentTimeMillis();
			for ( int i = 0 ;i < 10000; i++){
				double fourier = uFourier();
			}
			long t2_fourier = System.currentTimeMillis();
			System.out.println("delta_t_fourier: "+ (t2_fourier- t1_fourier));
		}
		
		double real = uReal();
		double fourier = uFourier();
		double self = uSelf();

		double nonReducedEnergy = real + fourier + self;
			
		if (false) { 
		System.out.println("total:               "+ nonReducedEnergy);
		System.out.println("real: "+ real);
		System.out.println("fourier: "+ fourier);
		System.out.println("self: "+ self);
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

	public void setBox(Box box) {
		// do nothing
		
	}
	//******************************** inner class ********************************************//
	public static class MyCharge{
		public MyCharge(double charge){
			this.charge = charge;
		}
		public final double charge;
		
	}
	
	
}
