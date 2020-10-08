/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.action.CalcVibrationalModes;
import etomica.action.WriteConfiguration;

import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.box.storage.Tokens;
import etomica.box.storage.VectorStorage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorBox;
import etomica.molecule.IMoleculeList;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

import java.io.FileWriter;
import java.io.IOException;


/**
 * 	Henkelman's Dimer Method (Rotation and Translation).
 * 
 * 	 | Finding saddle points on high dimensional surfaces.
 * 	 |    J. Chem. Phys., Vol. 111, No. 15, 15 October 1999.
 * 
 * 	@author msellers
 */

public class IntegratorDimerRT extends IntegratorBox implements AgentSource<Vector> {

	public final Box box1, box2;
	public Simulation sim;
	public double deltaR;
	public double dTheta, deltaXl, dXl;
	public double deltaTheta;
	public double curvature;
	public double deltaXmax;
	public double Frot, Frot2, Fprimerot, dFrot;
	public double saddleT;
	public double dFsq;
	public double gammai;
	public double saddleEnergy;
	public int rotNum;
	public int counter, rotCounter;
	public int movableAtoms;
	public boolean rotate, ortho, ortho2, startOrtho, saddleFound;
	public MeterPotentialEnergy energyBox1, energyBox2, energyBox0;
	public Vector[] THETA, THETAstar, THETAstarstar;
	public Vector[] F, F1, F2;
	public Vector[] Fperp, F1perp, F2perp;
	public Vector[] Gperp, Gperplast, Fperplast;
	public Vector[] Fstar, F1star, F2star, Fstarperp;
	public Vector[] Feff, Feffstar, Fr, Fpara;
	public Vector[] deltaV, V;
	public Vector[] newPosition;
	public Vector[] workVector3;
	public Vector[] N, Nstar, Neff, N1;
	public Vector NDelta, NstarDelta;
	public Vector workVector;
	public IRandom random1;
	public ISpecies [] movableSpecies;
	public PotentialCalculationForceSum force0, force1, force2;
	public AtomArrayList list, list1, list2;
	public final VectorStorage forces0, forces1, forces2;
	public IteratorDirective allatoms;
	public String file;
	public CalcVibrationalModes vib;

	
	public IntegratorDimerRT(Simulation sim, PotentialMaster potentialMaster,
							 ISpecies[] species, Box box) {
		this(sim, potentialMaster, sim.getRandom(), 1.0, species, box);
	}
	
	public IntegratorDimerRT(Simulation aSim, PotentialMaster potentialMaster,
							 IRandom random, double temperature,
							 ISpecies[] aspecies, Box box) {
		super(potentialMaster, temperature, box);
		this.random1 = random;
		this.sim = aSim;
		box1 = sim.makeBox(box.getBoundary());
		box2 = sim.makeBox(box.getBoundary());

		this.forces0 = box.getAtomStorage(Tokens.FORCES);
		this.forces1 = box1.getAtomStorage(Tokens.FORCES);
		this.forces2 = box2.getAtomStorage(Tokens.FORCES);
		this.force0 = new PotentialCalculationForceSum(forces0);
		this.force1 = new PotentialCalculationForceSum(forces1);
		this.force2 = new PotentialCalculationForceSum(forces2);
		this.allatoms = new IteratorDirective();
		this.movableSpecies = aspecies;

		deltaR = 1E-3;
		dXl = 1E-3;
		deltaXl = 0;
		deltaXmax = 0.04;
		
		deltaTheta = 0;
		dTheta = 1E-4;
		
		dFsq = 0.01;
		
		Frot = 1.0;
		dFrot = 0.1;
		
		counter = 0;
		rotCounter = 0;
		rotNum = 2;
		
		startOrtho = false;
		ortho2 = false;
	}
	
	public void setRotNum(int num){
		this.rotNum = num;
	}
	
	/**
	 * Set's the Dimer search's orthogonality criteria.
	 * @param orthoON (if true) Once the rotation scheme has converged to a small enough Frot value, the orthogonal search is started.
	 * @param o2 (if true) Launches another orthogonal path search after the primary orthogonal search.
	 * 
	 */
	public void setOrtho(boolean orthoON, boolean o2){
		startOrtho = orthoON;
		ortho2 = o2;
	}
	
	/**
     * Set's the filename of the output of this integrator.  "-minimum" is appended
     * to the string.
     * 
     * @param fileName String
     */
	public void setFileName(String fileName){
		file = fileName;
	}


	protected void doStepInternal() {
		
		//System.out.println("rotating...");
		
		rotateDimerNewton();
        
        //System.out.println("translating...");
        		
		translateDimerQuickmin();
		
		counter++;
	}
	
	public void reset() {
	    super.reset();
        
	    deltaR = 1E-3;
        dXl = 1E-3;
        deltaXl = 0;
        deltaXmax = 0.04;
        
        deltaTheta = 0;
        dTheta = 1E-4;
        
        dFsq = 0.05;
        
        Frot = 1.0;
        dFrot = 0.1;
        
        counter = 0;
        rotCounter = 0;
        rotNum = 1;

	    saddleFound=false;
        counter = 0;
        startOrtho = false;
        ortho2 = false;
	 // Use random unit array for N, N1 to generate dimer.

        // Normalize N, N1
        double mag = 0;
        double magN1 = 0;
        for (int i=0; i<N.length; i++){ 
            N1[i].setRandomSphere(random1);
            N[i].setRandomSphere(random1);
            mag += N[i].squared();
            magN1 += N1[i].squared();
        }
        mag = 1.0 / Math.sqrt(mag);
        magN1 = 1.0 / Math.sqrt(magN1);
        for (int i=0; i<N.length; i++){
            N[i].TE(mag);
            N1[i].TE(magN1);
        }

        // Set positions of atoms in replicas equal to box
        for(int i = 0; i<box.getLeafList().size(); i++){
            box1.getMoleculeList().get(i).getChildList().get(0).getPosition().E(box.getMoleculeList().get(i).getChildList().get(0).getPosition());
            box2.getMoleculeList().get(i).getChildList().get(0).getPosition().E(box.getMoleculeList().get(i).getChildList().get(0).getPosition());

        }
        // Offset replicas
        for(int i=0; i<N.length; i++){
            list1.get(i).getPosition().PEa1Tv1(deltaR, N[i]);
            list2.get(i).getPosition().PEa1Tv1(-deltaR, N[i]);
        }
        
        //System.out.println("...testing dimer direction.");
        if(energyBox0.getDataAsScalar()>energyBox1.getDataAsScalar()){
            //System.out.println(".S - Dimer pointed downhill, swapping ends.");
            for(int i=0; i<N.length; i++){
                  list1.get(i).getPosition().PEa1Tv1(-2.0*deltaR, N[i]);
                  list2.get(i).getPosition().PEa1Tv1(2.0*deltaR, N[i]);
              }
          }
        dimerForces(F1, F2, F);

	}
	
	/**
	 * Takes in a current configuration of atoms (Rc) and creates a dimer of their positions (R1 and R2).
	 * (R2)-------[Rc]-------(R1) ===> N
	 */
	protected void setup() {
		super.setup();

		movableAtoms = 0;
		for (int i = 0; i < movableSpecies.length; i++) {
			movableAtoms += box.getMoleculeList(movableSpecies[i]).size();
		}
		workVector = space.makeVector();

		N = new Vector[movableAtoms];
		Neff = new Vector[movableAtoms];
		Nstar = new Vector[movableAtoms];
		N1 = new Vector[movableAtoms];
		THETA = new Vector[movableAtoms];
		THETAstar = new Vector[movableAtoms];
		THETAstarstar = new Vector[movableAtoms];
		F = new Vector[movableAtoms];
		F1 = new Vector[movableAtoms];
		F2 = new Vector[movableAtoms];
		Fperp = new Vector[movableAtoms];
		Fperplast = new Vector[movableAtoms];
		Gperp = new Vector[movableAtoms];
		Gperplast = new Vector[movableAtoms];
		F1perp = new Vector[movableAtoms];
		F2perp = new Vector[movableAtoms];
		Fstar = new Vector[movableAtoms];
		F1star = new Vector[movableAtoms];
		F2star = new Vector[movableAtoms];
		Fstarperp = new Vector[movableAtoms];
		Feff = new Vector[movableAtoms];
		Feffstar = new Vector[movableAtoms];
		Fr = new Vector[movableAtoms];
		Fpara = new Vector[movableAtoms];
		deltaV = new Vector[movableAtoms];
		V = new Vector[movableAtoms];
		newPosition = new Vector[movableAtoms];
		workVector3 = new Vector[movableAtoms];

		for (int i = 0; i < movableAtoms; i++) {
			N[i] = space.makeVector();
			Neff[i] = space.makeVector();
			Nstar[i] = space.makeVector();
			N1[i] = space.makeVector();
			THETA[i] = space.makeVector();
			THETAstar[i] = space.makeVector();
			THETAstarstar[i] = space.makeVector();
			F[i] = space.makeVector();
			F1[i] = space.makeVector();
			F2[i] = space.makeVector();
			Fperp[i] = space.makeVector();
			Fperplast[i] = space.makeVector();
			Gperp[i] = space.makeVector();
			Gperplast[i] = space.makeVector();
			F1perp[i] = space.makeVector();
			F2perp[i] = space.makeVector();
			Fstar[i] = space.makeVector();
			F1star[i] = space.makeVector();
			F2star[i] = space.makeVector();
			Fstarperp[i] = space.makeVector();
			Feff[i] = space.makeVector();
			Feffstar[i] = space.makeVector();
			Fpara[i] = space.makeVector();
			deltaV[i] = space.makeVector();
			V[i] = space.makeVector();
			newPosition[i] = space.makeVector();
			workVector3[i] = space.makeVector();
		}

		if (potentialMaster instanceof PotentialMasterList) {
			this.getEventManager().addListener(((PotentialMasterList) potentialMaster).getNeighborManager(box1));
		}

		energyBox0 = new MeterPotentialEnergy(potentialMaster, box);
		energyBox1 = new MeterPotentialEnergy(potentialMaster, box1);
		energyBox2 = new MeterPotentialEnergy(potentialMaster, box2);

		for (int i = 0; i < sim.getSpeciesCount(); i++) {
			ISpecies species = sim.getSpecies(i);
			box1.setNMolecules(species, box.getNMolecules(species));
			box2.setNMolecules(species, box.getNMolecules(species));
		}

		// Atom list for movable and offset atoms
		list = new AtomArrayList();
		list1 = new AtomArrayList();
		list2 = new AtomArrayList();

		for (int i = 0; i < movableSpecies.length; i++) {
			IMoleculeList molecules = box.getMoleculeList(movableSpecies[i]);
			IMoleculeList molecules1 = box1.getMoleculeList(movableSpecies[i]);
			IMoleculeList molecules2 = box2.getMoleculeList(movableSpecies[i]);
			for (int j = 0; j < molecules.size(); j++) {
				list.add(molecules.get(j).getChildList().get(0));
				list1.add(molecules1.get(j).getChildList().get(0));
				list2.add(molecules2.get(j).getChildList().get(0));
			}
		}

	}
		

    /**
     * Rotates the dimer in potential energy hyperspace and aligns it with the surface's lowest 
     * curvature mode.  Using Newton's finite difference method, the rotational force perpendicular to the dimer's
     * orientation is minimized.  Minimization criteria is dependent on dFrot and rotCounter.
     * 
     * 
     */
	public void rotateDimerNewton(){
		rotCounter=0;
		
	    // Calculate new F
        dimerForces(F1, F2, F);
        
        dimerSaddleTolerance();
				
	    while(true){
			
	        // Calculate F1
            dimerForcesStar(F1, F2, F);
	        
	    	//ORTHOGONAL DIMER SEARCH
	    	if(ortho){
				//Remove N1 component from F's and orthogonal N
	    		// F1-(F1[dot]N1)*N1; F2-(F2[dot]N1)*N1; N-(N[dot]N1)*N1
	    		double magF1 = 0;
	    		double magF2 = 0;
	    		double magN1 = 0;
	    		for(int i=0; i<F1.length; i++){
	    	        magF1 += F1[i].dot(N1[i]);
	    	        magF2 += F2[i].dot(N1[i]);
	    			magN1 += N[i].dot(N1[i]);
	    		}
	    		for (int i=0; i<F1.length; i++){
	    			F1[i].PEa1Tv1(-magF1, N1[i]);	
	    			F2[i].PEa1Tv1(-magF2, N1[i]);
	    			N[i].PEa1Tv1(-magN1,N1[i]);
	    		}
        		System.out.println("ortho N: "+N[0].getX(0)+"    "+N[0].getX(1)+"    "+N[0].getX(2));
        		System.out.println(" lcm N1: "+N1[0].getX(0)+"    "+N1[0].getX(1)+"    "+N1[0].getX(2));
        		
        		//N is now a vector orthogonal to N1, check curvatures along each;
                if(ortho){
                    testOrthoCurvature(N1, N);
                }
			}
	    	
			// Calculate F|_
			dimerForcePerp(N, F1, F2, Fperp);
			
			/**
			//CONJUGATE GRADIENT SEARCH FOR ROTATIONAL PLANE
			// Find Gperp = Fperp + gammai*|Gperplast|*THETAstarstar
			//gammaI = ( (Fperpi - Fperpi-1)[dot]Fperpi )/(Fperpi[dot]Fperpi)
			
			double gt = 0.0;
			double gb = 0.0;
			for(int i=0; i<Fperp.length; i++){
				Fperplast[i].TE(-1.0);
				Fperplast[i].PE(Fperp[i]);
				gt += Fperplast[i].dot(Fperp[i]);
			}
			for(int i=0; i<Fperp.length; i++){
				gb += Fperp[i].dot(Fperp[i]);
			}
			gammai = gt/gb;
			double magG=0;
			for(int i=0; i<Gperplast.length; i++){
				magG += Gperplast[i].squared();
			}
			magG = Math.sqrt(magG);
			for(int i=0; i<Gperp.length; i++){
				Gperp[i].E(THETAstarstar[i]);
				Gperp[i].TE(gammai*magG);
				Gperp[i].PE(Fperp[i]);
			}			
			if(rotCounter==0){
				for(int i=0; i<Gperp.length; i++){
					Gperp[i].E(Fperp[i]);
				}
			}
			*/
			// Find THETA and normalize
			double mag = 0;
			for(int i=0; i<Fperp.length; i++){
				mag += Fperp[i].squared();
			}
			mag = 1.0 / Math.sqrt(mag);		
			for(int i=0; i<THETA.length; i++){
			    //THETA[i].E(Gperp[i]);
			    THETA[i].E(Fperp[i]);
				THETA[i].TE(mag);
			}
			//Copy arrays to remember for next iteration
			for(int i=0; i<Fperplast.length; i++){
				Gperplast[i].E(Gperp[i]);
				Fperplast[i].E(Fperp[i]);
				THETAstarstar[i].E(THETA[i]);
			}
			
            // NEWTON'S LINE MINIMIZATION FOR ROTATIONAL FORCE	
			Frot = 0;
            for(int i=0; i<Fperp.length; i++){
                  //Frot += Gperp[i].dot(THETA[i]);
                  Frot += Fperp[i].dot(THETA[i]);
            }
            
            // Leave loop if Frot at current location is small
            if(Frot<dFrot){
            	
            	//First time after converging to lowest curvature mode, N is copied to N1
            	//Start ortho
            	if(startOrtho){
            		System.out.println("Starting othogonal dimer search.");
            		
            		Vector workVec = space.makeVector();
            		for(int i=0;i<N.length;i++){
            			//swap
            			workVec.E(N1[i]);
            			N1[i].E(N[i]);
            			N[i].E(workVec);
            		}
            		
            		//Make N orthogonal to N1:  N - (N[dot]N1)*N1
            		double orthN = 0;
            		for(int i=0; i<N1.length; i++){
            			orthN += N[i].dot(N1[i]);
            		}
            		for (int i=0; i<N1.length; i++){
            			N[i].PEa1Tv1(-orthN,N1[i]);
            		}
            		System.out.println("ortho N: "+N[0].getX(0)+"    "+N[0].getX(1)+"    "+N[0].getX(2));
            		System.out.println(" lcm N1: "+N1[0].getX(0)+"    "+N1[0].getX(1)+"    "+N1[0].getX(2));
            		ortho=true;
            		startOrtho = false;
            		file = file+"_ortho";
            	}
            	
            	// Calculate curvature value
        		dimerCurvature(N, F1, F2);
        		
            	//System.out.println(rotCounter+" rotations.");
            	break;
            }
            
            double sinDtheta = Math.sin(dTheta);
            double cosDtheta = Math.cos(dTheta);
            
			// Find Nstar after dTheta rotation
	        Vector workVectorN1 = space.makeVector();
			for(int i=0; i<N.length; i++){
				workVectorN1.Ea1Tv1(cosDtheta, N[i]);
				workVectorN1.PEa1Tv1(sinDtheta, THETA[i]);
				Nstar[i].E(workVectorN1);
				workVectorN1.Ea1Tv1(-sinDtheta, N[i]);
				workVectorN1.PEa1Tv1(cosDtheta, THETA[i]);
				THETAstar[i].E(workVectorN1);
			}
			// Use N* to offset(rotate) replicas
			for(int i=0; i<Nstar.length; i++){
			    workVector.E(list.get(i).getPosition());
			    workVector.PEa1Tv1(deltaR, Nstar[i]);
                list1.get(i).getPosition().E(workVector);
                
                workVector.E(list.get(i).getPosition());
                workVector.PEa1Tv1(-deltaR, Nstar[i]);
                list2.get(i).getPosition().E(workVector);
            }
			// Calculate F*'s
			dimerForcesStar(F1star, F2star, F);     
			
			// Calculate F*|_
			dimerForcePerp(Nstar, F1star, F2star, Fstarperp);
						
			// second part of Frot
			for(int i=0; i<Fstarperp.length; i++){
				Frot += Fstarperp[i].dot(THETAstar[i]);
			}
			Frot = Frot/2.0;
			
			Fprimerot = 0;
			for(int i=0; i<Fstarperp.length; i++){
				Fprimerot += Fstarperp[i].dot(THETAstar[i]);
			}
			for(int i=0; i<Fperp.length; i++){
				//Fprimerot -= Gperp[i].dot(THETA[i]);
			    Fprimerot -= Fperp[i].dot(THETA[i]);
			}
			Fprimerot = Fprimerot / dTheta;
			
			// Find actual rotation angle to minimize energy
			deltaTheta = -0.5 * Math.atan(2.0*Frot/Fprimerot) - dTheta/2.0;				
			if(Fprimerot>0){deltaTheta = deltaTheta + Math.PI/2.0;}
			
			//System.out.println(".R - DeltaTheta "+deltaTheta+"  Frot "+Frot+"  Fprimerot "+Fprimerot);
            
			double sindeltaTheta = Math.sin(deltaTheta);
            double cosdeltaTheta = Math.cos(deltaTheta);
			
            /*
            if(deltaTheta < 0.1 && dTheta > 10E-7){
            	dTheta = dTheta/100.0;
            	System.out.println("R. - Scaling dTheta to "+dTheta);
            }
            */
			
            // Find N**
            Vector workVectorN2 = space.makeVector();
            for(int i=0; i<N.length; i++){               
                workVectorN2.Ea1Tv1(cosdeltaTheta, Nstar[i]);
                workVectorN2.PEa1Tv1(sindeltaTheta, THETAstar[i]);
                N[i].E(workVectorN2);
            }
            
            // Use new N to offset(rotate) replicas
            for(int i=0; i<N.length; i++){             
                workVector.E(list.get(i).getPosition());
                workVector.PEa1Tv1(deltaR, N[i]);
                list1.get(i).getPosition().E(workVector);
                
                workVector.E(list.get(i).getPosition());
                workVector.PEa1Tv1(-deltaR, N[i]);
                list2.get(i).getPosition().E(workVector);
            }        
            
			rotCounter++;
			
			if(rotCounter>rotNum){
				//System.out.println(".R - "+rotCounter+" rotations.");
				
				// Calculate estimated Curvature
				dimerCurvature(Nstar, F1star, F2star);
				double c1 = 0;
				for(int i=0; i<Fstarperp.length; i++){
					c1 += Fstarperp[i].squared();
				}
				c1 = 1.0 / Math.sqrt(c1);
				//System.out.println(" curveI "+curvature);
				curvature = curvature - (0.5 * c1 * Math.tan(deltaTheta - 0.5*dTheta));
				//System.out.println(" curveE "+curvature);
				break;
			}
		}
	}
	
	/**
	 * Moves the dimer along the lowest curvature mode of the energy surface, Quickmin style.
	 * Step length is determined using a line search (Newton) and the step direction is obtained
	 * by a calculation of the effective force on the dimer.
	 *  
	 */
	public void translateDimerQuickmin(){
			
		// Calculate Feff
		dimerForceEff(F, Feff, N);
		
		// Calculate Effective force normal
		double mag = 0;
		for(int i=0; i<Feff.length; i++){
		    mag += Feff[i].squared();
		    Neff[i].E(Feff[i]);
		}
		mag = 1.0/Math.sqrt(mag);
        for (int i=0; i<N.length; i++){
            Neff[i].TE(mag);
        }
        
        // Calculate magnitude before step
        double Feffmag = 0;
        for(int i=0; i<Feff.length; i++){
            Feffmag += Feff[i].dot(Neff[i]);
        }
        //System.out.println(".T - Feffmag: "+Feffmag);

        
		// Line search for new step length, move atoms a test distance
        dimerUpdatePositions(dXl, Neff);
		
		// Calculate new dimer center force
		dimerForceCenter(F);
        
		dimerSaddleTolerance();
		
		// Calculate Feffstar
        dimerForceEff(F, Feffstar, N);
		
        // Calculate magnitude of step
        Feffmag = 0;
        for(int i=0; i<Feffstar.length; i++){
            workVector3[i].Ev1Pv2(Feffstar[i],Feff[i]);
            Feffmag += workVector3[i].dot(Neff[i]);
        }
        Feffmag = Feffmag/2.0;
        
        // Calculate new curvature
        double curvatureLM = 0;
        for(int i=0; i<Feffstar.length; i++){
            workVector3[i].Ev1Mv2(Feffstar[i],Feff[i]);
        }
        for(int i=0; i<Feffstar.length; i++){
            curvatureLM += workVector3[i].dot(Neff[i]);
        }
        curvatureLM = curvatureLM/dXl;
        
        //System.out.println(".T - Feffmag: "+Feffmag+"  CurvatureLM: "+curvatureLM);
               
        deltaXl = ((-Feffmag/curvatureLM) + dXl/2.0);
        if(deltaXl<0.0 && Math.abs(deltaXl)>deltaXmax){deltaXl = -deltaXmax + dXl/2.0;}
        if(deltaXl>deltaXmax || curvature>0.0){
            deltaXl = deltaXmax-(dXl/2.0);
        }
        
        //System.out.println(".T - Step calculated at "+deltaXl+".  curvature is "+curvature);
        
        // Step dimer
        dimerUpdatePositions(deltaXl, Neff);
              
        // Calculate new Normal
        dimerNormal();    
	}
	
	/**
	 * Computes curvature value, C, for the energy surface.  Quantifies the difference
	 * in forces at each end of the dimer into a single number.
	 * 
	 */
	protected void dimerCurvature(Vector[] aN, Vector[] aF1, Vector[] aF2){
		curvature = 0.0;
		
	      // Copy forces of dimer end and center (R1, R) to local array
        for(int i=0; i<aF1.length; i++){
            aF2[i].Ea1Tv1(2.0, F[i]);
            aF2[i].ME(aF1[i]);  
        }
		
		// (F2 - F1)[dot]N / 2*deltaR
		for(int i=0; i<aF1.length; i++){
			workVector.E(aF2[i]);
			workVector.ME(aF1[i]);
			
			curvature += workVector.dot(aN[i]);
		}
		
		curvature = 0.5 * curvature / deltaR;
		
		//System.out.println(curvature+" curvature");
	}
	
	/**
	 * Checks the curvature value at dimer orientation along a vector orthogonal to
	 * lowest curvature mode's (lcm) vector with the lcm's unit vector.  When 
	 * the ortho vector returns a curvature less than that of the lcm vector, 
	 * ORTHO is turned off.
	 * 
	 */
	protected void testOrthoCurvature(Vector[] aN1, Vector[] aNortho){
		
		//Compute current curvature along orthogonal N
		dimerCurvature(aNortho, F1, F2);
		double orthoCurve = curvature;
		
		dimerCurvature(aN1, F1, F2);
		double n1Curve = curvature;
		
		System.out.println("curvatures: "+n1Curve+" (LCM)    "+orthoCurve+" (ORTHO)");
		
		if(orthoCurve<n1Curve){
			System.out.println("Secondary mode found, stopping ortho search.");
			ortho = false;
			
			// Run another Orthogonal search (N1, N2)
			if(ortho2){
				startOrtho = true;
				ortho2=false;
			}
		}

	}
	
	// Compute Normal vector for dimer orientation
	protected void dimerNormal(){
	    double mag=0;
		Vector workvector;
		workvector = space.makeVector();
		
		// N =  (R1 - R2) / (-2*deltaR)
		for (int i=0; i<N.length; i++){	
			workvector.E(list1.get(i).getPosition());
			workvector.ME(list2.get(i).getPosition());
			N[i].E(workvector);
			mag += workvector.squared();
		}
		mag = 1.0/Math.sqrt(mag);		
		for (int i=0; i<N.length; i++){
		    N[i].TE(mag);
		}
	}
		
	// Reset forces in boxes 1 and 2, call calculate, and copy over new forces
	protected void dimerForces(Vector[] aF1, Vector[] aF2, Vector[] aF){
		force1.reset();
		force0.reset();
		
		potentialMaster.calculate(box1, allatoms, force1);
		potentialMaster.calculate(box, allatoms, force0);
		
		// Copy forces of dimer end and center (R1, R) to local array
		for(int i=0; i<aF1.length; i++){
			aF1[i].E(forces1.get(list1.get(i).getLeafIndex()));
			aF[i].E(forces0.get(list.get(i).getLeafIndex()));
			aF2[i].Ea1Tv1(2.0, aF[i]);
			aF2[i].ME(aF1[i]);	
		}
	}
	
	
	protected void dimerForcesStar(Vector[] aF1star, Vector[] aF2star, Vector[] aF){
	    force1.reset();
	    potentialMaster.calculate(box1, allatoms, force1);
	    
	 // Copy forces of dimer end and center (R1, R) to local array
	    for(int i=0; i<aF1star.length; i++){
			aF1star[i].E(forces1.get(list1.get(i).getLeafIndex()));
			aF2star[i].Ea1Tv1(2.0, aF[i]);
			aF2star[i].ME(aF1star[i]);	
		}
	}
	
	protected void dimerForceCenter(Vector[] aF){
	    force0.reset();
	    potentialMaster.calculate(box, allatoms, force0);
	    
	 // Copy forces of dimer end and center (R1, R) to local array
	    for(int i=0; i<aF.length; i++){
	    	aF[i].E(forces0.get(list.get(i).getLeafIndex()));
		}
	}
	
	
	// Calculate the force, aF, perpendicular to dimer orientation N.
	protected void dimerForcePerp(Vector[] aN, Vector[] aF1, Vector[] aF2, Vector[] aFperp){
		
		double mag1 = 0;
		double mag2 = 0;
		
		// F|_ = F1 - (F1[dot]N)*N  - F2  +  (F2[dot]N)*N
		for(int i=0; i<aF1.length; i++){
			mag1 += aF1[i].dot(aN[i]);
			mag2 += aF2[i].dot(aN[i]);
		}
		for (int i=0; i<aFperp.length; i++){
			aFperp[i].E(aF1[i]);
			aFperp[i].PEa1Tv1(-mag1, aN[i]);	
			aFperp[i].ME(aF2[i]);
			aFperp[i].PEa1Tv1(mag2, aN[i]);
			aFperp[i].TE(1.0/deltaR);
		}

	}
		
	// Calculate the effective force on the dimer using curvature and F, N
	protected void dimerForceEff(Vector[] aF, Vector[] aFeff, Vector[] aNeff){
		
	    double mag = 0;
	    for(int i=0; i<aF.length; i++){
            mag += aF[i].dot(aNeff[i]);
        }
	    //System.out.println(F[0]+" actual force");
	    //System.out.println(mag+" F dot N");
	    
		// Feff = F - 2(F[dot]N)N
		if(curvature<0.0){
			Vector workvector = space.makeVector();
			for (int i=0; i<aNeff.length; i++){
			    workvector.Ea1Tv1(2.0*mag, aNeff[i]);
			    aFeff[i].Ev1Mv2(aF[i], workvector);
			}
		}
		
		// Feff = -(F[dot]N)N
		else{
			for(int i=0; i<aFeff.length; i++){
				aFeff[i].Ea1Tv1(-mag, aNeff[i]);
			}
		}
	}
	
	// Update positions according Henkelman 2004
	protected void dimerUpdatePositions(double a1, Vector[] normal){
		Vector workvector = space.makeVector();
		//System.out.println(".T - Stepping "+a1);
		for(int i=0; i<normal.length; i++){
		    workvector.Ea1Tv1(a1, normal[i]);			
    		list.get(i).getPosition().PE(workvector);
    		list1.get(i).getPosition().PE(workvector);
    		list2.get(i).getPosition().PE(workvector);
    	}	

	}
	
	protected void resizeDimer(){
	    
	    dimerNormal();
	    
	    // Set positions of atoms in replicas equal to box
        for(int i = 0; i<box.getLeafList().size(); i++){
            box1.getMoleculeList().get(i).getChildList().get(0).getPosition().E(box.getMoleculeList().get(i).getChildList().get(0).getPosition());
            box2.getMoleculeList().get(i).getChildList().get(0).getPosition().E(box.getMoleculeList().get(i).getChildList().get(0).getPosition());

        }
        // Offset replicas
        for(int i=0; i<N.length; i++){
            list1.get(i).getPosition().PEa1Tv1(deltaR, N[i]);
            list2.get(i).getPosition().PEa1Tv1(-deltaR, N[i]);
        }
	}
	
	// Calculates and checks magnitude of 3N dimensional force vector
	protected void dimerSaddleTolerance(){
		saddleT = 0.0;
		//If every force is less than dF, consider saddle point found.
		for(int i=0; i<F.length; i++){
			saddleT += F[i].squared();
		}
		if(saddleT<0.5){
		    setRotNum(5);
		    if(energyBox0.getDataAsScalar()>energyBox1.getDataAsScalar()){
		        //deltaR = 1E-5;
                //resizeDimer();
		        dXl = 0.000001;
    		    deltaXmax = 0.001;
    		    dTheta = 0.000001;
    		    dFrot = 0.01;
		    }
		}
		if(saddleT<dFsq){
		    
		    vib = new CalcVibrationalModes();
	        vib.setup(box, potentialMaster, box.getMoleculeList(movableSpecies[0]), space);
	        vib.actionPerformed();
	        
	        
		    saddleFound = true;
		    saddleEnergy = energyBox0.getDataAsScalar();
		    
		    
	        try{
	            FileWriter fileWriter = new FileWriter(file+"_s_ev", false);
	            fileWriter.write(energyBox0.getDataAsScalar()+"\n"+vib.getProductOfFrequencies());
	            fileWriter.close();
	        }catch(IOException e) {
	          
	        }
	        		
			
			// Write out configurations of 3 boxes
		    WriteConfiguration writer = new WriteConfiguration(space);
		    writer.setConfName(file+"_saddle");
		    writer.setBox(box);
		    writer.actionPerformed();
		    
		    writer.setConfName(file+"_A_saddle");
		    writer.setBox(box1);
		    writer.actionPerformed();
		    
		    writer.setConfName(file+"_B_saddle");
		    writer.setBox(box2);
		    writer.actionPerformed();
		    
		    //activityIntegrate.setMaxSteps(0);
		    
		}
			
	}

	public Vector makeAgent(IAtom a, Box agentBox) {
		return space.makeVector();
	}

	public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {}
	
	/*
	public static class PotentialMasterListDimer extends PotentialMasterList{

		public PotentialMasterListDimer(Simulation sim, ISpace space) {
			super(sim, space);
			
		}
		
		public void setSpecies(ISpecies[] species){
			this.species = species;
		}
		
	   public void calculate(Box box, IteratorDirective id, PotentialCalculation pc) {
	        if(!enabled) return;
	        IAtom targetAtom = id.getTargetAtom();
	        if (targetAtom != null) {
	        	super.calculate(box, id, pc);
	        	return;
	        }
	        NeighborListManager neighborManager = (NeighborListManager)neighborListAgentManager.getAgent(box);
	        
            if (Debug.ON && id.direction() != IteratorDirective.Direction.UP) {
                throw new IllegalArgumentException("When there is no target, iterator directive must be up");
            }
            // invoke setBox on all potentials
            for (int i=0; i<allPotentials.length; i++) {
                allPotentials[i].setBox(box);
            }

            //no target atoms specified
            //call calculate with each SpeciesAgent
	        for(int j=0; j<species.length; j++){    
            	IAtomSet list = box.getMoleculeList(species[j]);
	            int size = list.getAtomCount();
	            for (int i=0; i<size; i++) {
	                calculate((IMolecule)list.getAtom(i), id, pc, neighborManager);//call calculate with the SpeciesAgent
	            }
	        }
	   }
		
		
		public ISpecies [] species;
	}
	*/
}
