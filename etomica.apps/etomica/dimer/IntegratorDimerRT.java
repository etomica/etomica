package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.species.Species;
import etomica.util.IRandom;

/**
 * 	Henkelman's Dimer Method (Rotation and Translation).
 * 
 * 	 | Finding saddle points on high dimensional surfaces.
 * 	 |    J. Chem. Phys., Vol. 111, No. 15, 15 October 1999.
 * 
 * 	@author msellers
 */

public class IntegratorDimerRT extends IntegratorBox implements AgentSource {

	public Box box1, box2;
	public ISimulation sim;
	public double deltaT;
	public double deltaR;
	public double dTheta, deltaXl, dXl;
	public double deltaTheta;
	public double curvature;
	public double deltaXmax;
	public double Frot, Frot2, Fprimerot, dFrot;
	public double sinDtheta, cosDtheta;
	public double saddleT;
	public double dFsq;
	public int counter, rotCounter;
	public int movableAtoms;
	public boolean rotate;
	public MeterPotentialEnergy energyBox1, energyBox2, energyBox0;
	public IVector [] THETA, THETAstar;
	public IVector [] F, F1, F2;
	public IVector [] Fperp, F1perp, F2perp;
	public IVector [] Fstar, F1star, F2star, Fstarperp;
	public IVector [] Feff, Feffstar, Fr, Fpara;
	public IVector [] deltaV, V;
	public IVector [] newPosition;
	public IVector [] workVector3;
	public IVectorRandom [] N, Neff;
	public IVector NDelta, NstarDelta;
	public IVector workVector1;
	public IVector workVector2;
	public IRandom random1;
	public Species [] movableSpecies;
	public PotentialCalculationForceSum force0, force1, force2;
	public AtomArrayList list, list1, list2;
	public AtomAgentManager atomAgent0, atomAgent1, atomAgent2;
	public IteratorDirective allatoms;
	public String file;
	public FileWriter fileWriter;
	public ActivityIntegrate activityIntegrate;
	
	
	public IntegratorDimerRT(ISimulation sim, PotentialMaster potentialMaster, Species[] species, String file) {
		this(sim, potentialMaster, sim.getRandom(), 0.005, 1.0, species, file);
	}
	
	public IntegratorDimerRT(ISimulation aSim, PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Species[] aspecies, String aFile) {
		super(potentialMaster, temperature);
		this.random1 = random;
		this.sim = aSim;
		this.force0 = new PotentialCalculationForceSum();
		this.force1 = new PotentialCalculationForceSum();
		this.force2 = new PotentialCalculationForceSum();
		this.allatoms = new IteratorDirective();
		this.deltaT = timeStep;
		this.movableSpecies = aspecies;
		this.file = aFile;
		
		
		deltaR = 10E-4;
		
		dXl = 10E-4;
		deltaXl = 0;
		deltaXmax = 0.05;
		dTheta = 10E-4;
		
		deltaTheta = 0;
		
		dFsq = 0.001*0.001;
		
		Frot = 1.0;
		dFrot = 0.1;
		
		counter = 0;
		rotCounter = 0;
		rotate = true;		
	}
	
	
	public void doStepInternal(){
	    	    
		rotateDimerNewton();
		rotCounter=0;
		translateDimerQuickmin();
		System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2));
		dimerSaddleTolerance();		
		counter++;
	}
			
	// Takes in a current configuration of atoms (Rc) and creates a dimer of their positions (R1 and R2).
	// (R2)-------[Rc]-------(R1) ===> N
	protected void setup() throws ConfigurationOverlapException{
		super.setup();
		
		movableAtoms = 0;
		for(int i=0; i<movableSpecies.length; i++){
		    movableAtoms += box.getMoleculeList(movableSpecies[i]).getAtomCount();
		}
		workVector1 = box.getSpace().makeVector();
        workVector2 = box.getSpace().makeVector();
		
        N = new IVectorRandom [movableAtoms];
        Neff = new IVectorRandom [movableAtoms];
        THETA = new IVector [movableAtoms];
        THETAstar = new IVector [movableAtoms];
        F = new IVector [movableAtoms];
        F1 = new IVector [movableAtoms];
        F2 = new IVector [movableAtoms];
        Fperp = new IVector [movableAtoms];
        F1perp = new IVector [movableAtoms];
        F2perp = new IVector [movableAtoms];
        Fstar = new IVector [movableAtoms];
        F1star = new IVector [movableAtoms];
        F2star = new IVector [movableAtoms];
        Fstarperp = new IVector [movableAtoms];
        Feff = new IVector [movableAtoms];
        Feffstar = new IVector [movableAtoms];
        Fr = new IVector [movableAtoms];
        Fpara = new IVector [movableAtoms];
        deltaV = new IVector [movableAtoms];
        V = new IVector [movableAtoms];
        newPosition = new IVector [movableAtoms];
        workVector3 = new IVector [movableAtoms];
        
        for(int i=0; i<movableAtoms; i++){
            N[i] = (IVectorRandom)box.getSpace().makeVector();
            Neff[i] = (IVectorRandom)box.getSpace().makeVector();
            THETA[i] = box.getSpace().makeVector();
            THETAstar[i] = box.getSpace().makeVector();
            F[i] = box.getSpace().makeVector();
            F1[i] = box.getSpace().makeVector();
            F2[i] = box.getSpace().makeVector();
            Fperp[i] = box.getSpace().makeVector();
            F1perp[i] = box.getSpace().makeVector();
            F2perp[i] = box.getSpace().makeVector();
            Fstar[i] = box.getSpace().makeVector();
            F1star[i] = box.getSpace().makeVector();
            F2star[i] = box.getSpace().makeVector();
            Fstarperp[i] = box.getSpace().makeVector();
            Feff[i] = box.getSpace().makeVector();
            Feffstar[i] = box.getSpace().makeVector();
            Fpara[i] = box.getSpace().makeVector();
            deltaV[i] = box.getSpace().makeVector();
            V[i] = box.getSpace().makeVector();
            newPosition[i] = box.getSpace().makeVector();
            workVector3[i] = box.getSpace().makeVector();
        }
		
		// Use random unit array for N to generate dimer.

		// Normalize N
		double mag = 0;
		for (int i=0; i<N.length; i++){ 
			N[i].setRandomSphere(random1);
			mag += N[i].squared();
		}
		
		mag = 1.0 / Math.sqrt(mag);
		for (int i=0; i<N.length; i++){
			N[i].TE(mag);
		}
		
		//Offset R1 and R2 (dimer ends) from initial configuration, along N.
		box1 = new Box(box.getBoundary());
		box2 = new Box(box.getBoundary());
		
		sim.addBox(box1);
		sim.addBox(box2);
		
		energyBox1 = new MeterPotentialEnergy(this.potential);
		energyBox2 = new MeterPotentialEnergy(this.potential);
		energyBox0 = new MeterPotentialEnergy(this.potential);
		energyBox1.setBox(box1);
		energyBox2.setBox(box2);
		energyBox0.setBox(box);
		
		atomAgent0 = new AtomAgentManager(this, box);
		atomAgent1 = new AtomAgentManager(this, box1);
		atomAgent2 = new AtomAgentManager(this, box2);
		
		force0.setAgentManager(atomAgent0);
		force1.setAgentManager(atomAgent1);
		force2.setAgentManager(atomAgent2);
		
		Species [] species = sim.getSpeciesManager().getSpecies();
		
		for(int i=0; i<species.length; i++){
			box1.setNMolecules(species[i], box.getNMolecules(species[i]));
			box2.setNMolecules(species[i], box.getNMolecules(species[i]));
		}
		
		// Set positions of atoms in replicas equal to box
		for(int i=0; i<box.atomCount(); i++){
			((IAtomPositioned)box1.getLeafList().getAtom(i)).getPosition().E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
			((IAtomPositioned)box2.getLeafList().getAtom(i)).getPosition().E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
		}
		
		// Atom list for movable and offset atoms
		list = new AtomArrayList();
		list1 = new AtomArrayList();
		list2 = new AtomArrayList();
		
		for(int i=0; i<movableSpecies.length; i++){
            AtomSet molecules = box.getMoleculeList(movableSpecies[i]);
            AtomSet molecules1 = box1.getMoleculeList(movableSpecies[i]);
            AtomSet molecules2 = box2.getMoleculeList(movableSpecies[i]);
            for (int j=0; j<molecules.getAtomCount(); j++) {
            	
                list.add(((IMolecule)molecules.getAtom(j)).getChildList().getAtom(0));
                list1.add(((IMolecule)molecules1.getAtom(j)).getChildList().getAtom(0));
                list2.add(((IMolecule)molecules2.getAtom(j)).getChildList().getAtom(0));
            }
		}
		
		// Offset replicas
		for(int i=0; i<N.length; i++){
			((IAtomPositioned)list1.getAtom(i)).getPosition().PEa1Tv1(deltaR, N[i]);
			((IAtomPositioned)list2.getAtom(i)).getPosition().PEa1Tv1(-deltaR, N[i]);
		}
		
		
		// Write out initial configuration
		System.out.println(file+" ***Dimer Saddle Search***");
		System.out.println(N[0]);
		System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2)
				+"     "+((IAtomPositioned)list1.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list1.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list1.getAtom(0)).getPosition().x(2)
				+"     "+((IAtomPositioned)list2.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list2.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list2.getAtom(0)).getPosition().x(2));
		
		// Calculate F's
		dimerForces(F1, F2, F);
	}
		
	// Rotate the dimer to align it with the lowest curvature mode of the potential energy surface, NEWTON STYLE.
	public void rotateDimerNewton(){
	    dTheta = 10E-5;
	    deltaTheta = 1.0;
	    Frot = 1.0;
	    
	    while(true){
			
			// Calculate F|_
			dimerForcePerp(N, F1, F2, Fperp);
			
			// Find THETA and normalize
			double mag = 0;
			for(int i=0; i<Fperp.length; i++){
				mag += Fperp[i].squared();
			}
			mag = 1.0 / Math.sqrt(mag);		
			    
			for(int i=0; i<THETA.length; i++){
				THETA[i].E(Fperp[i]);
				THETA[i].TE(mag);
			}
						
            // Compute scalar rotational force and change in rotational force over length dTheta
			Frot = 0;
            for(int i=0; i<Fperp.length; i++){
                  Frot += Fperp[i].dot(THETA[i]);
            }
            
           // Leave loop if Frot at current location is small
            if(Frot<dFrot){break;}
            
            sinDtheta = Math.sin(dTheta);
	        cosDtheta = Math.cos(dTheta);
			
	     // Find THETA* in rotational plane
            IVector workVectorTHETA1 = box.getSpace().makeVector();
            for(int i=0; i<N.length; i++){
                workVectorTHETA1.Ea1Tv1(cosDtheta, THETA[i]);
                workVectorTHETA1.PEa1Tv1(-sinDtheta, N[i]);
                THETAstar[i].E(workVectorTHETA1);
            }
	        
			// Find new N after dTheta rotation
	        IVector workVectorN1 = box.getSpace().makeVector();
			for(int i=0; i<N.length; i++){
				workVectorN1.Ea1Tv1(cosDtheta, N[i]);
				workVectorN1.PEa1Tv1(sinDtheta, THETA[i]);
				N[i].E(workVectorN1);
			}
						
			// Use new N to offset(rotate) replicas
			IVector workVector1 = box.getSpace().makeVector();
			for(int i=0; i<N.length; i++){
			    workVector1.E(((IAtomPositioned)list.getAtom(i)).getPosition());
			    workVector1.PEa1Tv1(deltaR, N[i]);
                ((IAtomPositioned)list1.getAtom(i)).getPosition().E(workVector1);
                
                workVector1.E(((IAtomPositioned)list.getAtom(i)).getPosition());
                workVector1.PEa1Tv1(-deltaR, N[i]);
                ((IAtomPositioned)list2.getAtom(i)).getPosition().E(workVector1);
            }
			
			// Calculate F*'s
			dimerForcesStar(F, F1star, F2star, Fstar);     
			
			// Calculate F*|_
			dimerForcePerp(N, F1star, F2star, Fstarperp);
			
			// second part of Frot
			for(int i=0; i<Fstarperp.length; i++){
				Frot += Fstarperp[i].dot(THETAstar[i]);
			}
			
			Frot = Frot/2.0;
			
			Fprimerot = 0;
			for(int i=0; i<Fstar.length; i++){
				Fprimerot += Fstarperp[i].dot(THETAstar[i]);
			}
			for(int i=0; i<Fstar.length; i++){
				Fprimerot -= Fperp[i].dot(THETA[i]);
			}
			
			Fprimerot = Fprimerot / dTheta;
			
			// Find actual rotation angle to minimize energy
			
			// Adatom on surface 
			deltaTheta = -0.5 * Math.atan(2.0*Frot/Fprimerot) - dTheta/2.0;
						
			if(Fprimerot>0){deltaTheta = deltaTheta + Math.PI/2.0;}
			
			//double energy1 = (energyBox1.getDataAsScalar()-energyBox0.getDataAsScalar());
			//double energy2 = (energyBox2.getDataAsScalar()-energyBox0.getDataAsScalar());
			// System.out.println(energy1+energy2+"     "+Frot+"     "+Fprimerot+"     "+deltaTheta+"     "+N[0]);
			
			 // Check deltaTheta vs. dTheta and adjust step size
			
			if (deltaTheta < 0){
                    dTheta /= 100;
            //        System.out.println("dTheta changed to "+dTheta);
                    if (deltaTheta < -dTheta){
                        deltaTheta = -dTheta;
                    }
            }
						
			double sindeltaTheta = Math.sin(deltaTheta);
            double cosdeltaTheta = Math.cos(deltaTheta);
			
			// Find N**
            IVector workVectorN2 = box.getSpace().makeVector();
            for(int i=0; i<N.length; i++){               
                workVectorN2.Ea1Tv1(cosdeltaTheta, N[i]);
                workVectorN2.PEa1Tv1(sindeltaTheta, THETAstar[i]);
                N[i].E(workVectorN2);

            }
            
            // Use new N to offset(rotate) replicas
            IVector workVector2 = box.getSpace().makeVector();
            for(int i=0; i<N.length; i++){             
                workVector2.E(((IAtomPositioned)list.getAtom(i)).getPosition());
                workVector2.PEa1Tv1(deltaR, N[i]);
                ((IAtomPositioned)list1.getAtom(i)).getPosition().E(workVector2);
                
                workVector2.E(((IAtomPositioned)list.getAtom(i)).getPosition());
                workVector2.PEa1Tv1(-deltaR, N[i]);
                ((IAtomPositioned)list2.getAtom(i)).getPosition().E(workVector2);
            }        
            
			// Calculate new F's
			dimerForces(F1, F2, F);
			
			// Calculate new Normal
			dimerNormal();
						    
			rotCounter++;
			
			if(rotCounter>15){
			    break;
			}
		}
	    
	}
	
	// Moves the dimer along the lowest curvature mode of the energy surface, Quickmin style.  Runs after rotateDimer[X]().
	public void translateDimerQuickmin(){
	
		// Calculate curvature value
		dimerCurvature(N, F1, F2);
			
		// Calculate Feff
		dimerForceEff(F, Feff);
		
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
        
		// Line search for new step length, move atoms a test distance
		for(int i=0; i<N.length; i++){
		    newPosition[i].Ea1Tv1(dXl,Neff[i]);
		    ((IAtomPositioned)list.getAtom(i)).getPosition().PE(newPosition[i]);
		    ((IAtomPositioned)list1.getAtom(i)).getPosition().PE(newPosition[i]);
		    ((IAtomPositioned)list2.getAtom(i)).getPosition().PE(newPosition[i]);
		}
		
		// Calculate new Fstar
		dimerForces(F1, F2, F);
		
		// Calculate new Normal
        dimerNormal();
        		
		// Calculate Feffstar
        dimerForceEff(F, Feffstar);
		
        // Calculate magnitude of step
        double Feffmag = 0;
        for(int i=0; i<Feffstar.length; i++){
            workVector3[i].Ev1Pv2(Feffstar[i],Feff[i]);
        }
        for(int i=0; i<Feffstar.length; i++){
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
               
        deltaXl = ((-Feffmag/curvatureLM) - (dXl/2));
        
        //System.out.println(Feff[0]+" effective force");
        //System.out.println(Neff[0]+" effective normal");
        //System.out.println(N[0]+" dimer normal");
        
        // Step dimer
        dimerUpdatePositions(deltaXl, Neff);
        
        // Calculate new F
        dimerForces(F1, F2, F);
        
        // Calculate new Normal
        dimerNormal();
        
	}
	
	// Compute curvature value, C, for the energy surface
	protected void dimerCurvature(IVectorRandom [] aN, IVector [] aF1, IVector [] aF2){
		curvature = 0.0;
		
		// (F2 - F1)[dot]N / 2*deltaR
		for(int i=0; i<aF1.length; i++){
			workVector1.E(aF2[i]);
			workVector1.ME(aF1[i]);
			
			curvature += workVector1.dot(aN[i]);
		}
		
		curvature = 0.5 * curvature / deltaR;
		
		//System.out.println(curvature+" curvature");
	}
	
	// Compute Normal vector for dimer orientation
	protected void dimerNormal(){
	    double mag=0;
		IVector workvector;
		workvector = box.getSpace().makeVector();
		
		// N =  (R1 - R2) / (-2*deltaR)
		for (int i=0; i<N.length; i++){	
			workvector.E(((IAtomPositioned)list1.getAtom(i)).getPosition());
			workvector.ME(((IAtomPositioned)list2.getAtom(i)).getPosition());
			N[i].E(workvector);
			mag += workvector.squared();
		}
		
		mag = 1.0/Math.sqrt(mag);
		
		for (int i=0; i<N.length; i++){
		    N[i].TE(mag);
		}
	}
		
	// Reset forces in boxes 1 and 2, call calculate, and copy over new forces
	protected void dimerForces(IVector [] aF1, IVector [] aF2, IVector [] aF){
		force1.reset();
		force0.reset();
		
		potential.calculate(box1, allatoms, force1);
		potential.calculate(box, allatoms, force0);
		
		// Copy forces of dimer end and center (R1, R) to local array
		for(int i=0; i<aF1.length; i++){
			aF1[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent1.getAgent(list1.getAtom(i))).force());
			aF[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent0.getAgent(list.getAtom(i))).force());
			aF2[i].Ea1Tv1(2.0, aF[i]);
			aF2[i].ME(aF1[i]);
			
		}
	}
	
	
	protected void dimerForcesStar(IVector [] aF, IVector [] aF1star, IVector [] aF2star, IVector []aFstar){
	    force1.reset();
	    potential.calculate(box1, allatoms, force1);
	    
	 // Copy forces of dimer end and center (R1, R) to local array
        for(int i=0; i<aF1star.length; i++){
            aF1star[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent1.getAgent(list1.getAtom(i))).force());
            aFstar[i].E(aF[i]);
            aF2star[i].Ea1Tv1(2.0, aFstar[i]);
            aF2star[i].ME(aF1star[i]);
        }
	}
	
	// Calculate the force, aF, perpendicular to dimer orientation N.
	protected void dimerForcePerp(IVectorRandom [] aN, IVector [] aF1, IVector [] aF2, IVector [] aFperp){
		
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
	protected void dimerForceEff(IVector [] aF, IVector [] aFeff){
		
	    double mag = 0;
	    for(int i=0; i<aF.length; i++){
            mag += aF[i].dot(N[i]);
        }
	    //System.out.println(F[0]+" actual force");
	    //System.out.println(mag+" F dot N");
	    
		// Feff = F - 2(F[dot]N)N
		if(curvature<0){
			IVector workvector = box.getSpace().makeVector();
			for (int i=0; i<N.length; i++){
			    workvector.Ea1Tv1(2.0*mag, N[i]);
			    aFeff[i].Ev1Mv2(aF[i], workvector);
			}
		}
		
		// Feff = -(F[dot]N)N
		else{
			for(int i=0; i<aFeff.length; i++){
				aFeff[i].Ea1Tv1(-mag, N[i]);
			}
		}
	}
	
	// Update positions according Henkelman 2004
	protected void dimerUpdatePositions(double a1, IVector [] normal){
		
		IVector workvector = box.getSpace().makeVector();
		double sum = 0;
		
		if(a1>deltaXmax||curvature>0){
            a1 = deltaXmax - dXl;
        }
		workvector.Ea1Tv1(a1, normal[0]);
		//System.out.println(workvector+" step vector");
		
		for(int i=0; i<normal.length; i++){
		    workvector.Ea1Tv1(a1, normal[i]);
			
    		((IAtomPositioned)list.getAtom(i)).getPosition().PE(workvector);
    		((IAtomPositioned)list1.getAtom(i)).getPosition().PE(workvector);
    		((IAtomPositioned)list2.getAtom(i)).getPosition().PE(workvector);
    	}	

	}
	
	// Calculates and checks magnitude of 3N dimensional force vector
	protected void dimerSaddleTolerance(){
		saddleT = 0.0;
		//If every force is less than dF, consider saddle point found.
		for(int i=0; i<F.length; i++){
			saddleT += F[i].squared();
		}
		saddleT /= Math.sqrt(saddleT);
		if(saddleT<dFsq){
			System.out.println(file+" +++Dimer Saddle Found+++");
			System.out.println("  -"+counter+" steps.  "+saddleT+" magnitude of force array.  "+energyBox0.getDataAsScalar()+" energy of box 0.");

			// Write out configurations of 3 boxes
		    WriteConfiguration writer = new WriteConfiguration();
		    writer.setConfName(file+"_saddle");
		    writer.setBox(box);
		    writer.actionPerformed();
		    
		    writer.setConfName(file+"_A_saddle");
		    writer.setBox(box1);
		    writer.actionPerformed();
		    
		    writer.setConfName(file+"_B_saddle");
		    writer.setBox(box2);
		    writer.actionPerformed();
		    
	        activityIntegrate.setMaxSteps(0);
			//System.exit(1);
		}
			
	}
	
	public void setActivityIntegrate(ActivityIntegrate ai){
		activityIntegrate = ai;
	}
	
	public Class getAgentClass() {
		return IntegratorVelocityVerlet.MyAgent.class;
	}

	public Object makeAgent(IAtom a) {
		return new IntegratorVelocityVerlet.MyAgent(box.getSpace());
	}

	public void releaseAgent(Object agent, IAtom atom) {
		// TODO Auto-generated method stub	
	}
	

}
