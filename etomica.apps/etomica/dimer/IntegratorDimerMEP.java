package etomica.dimer;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
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
 * 	Henkelman's Dimer Method (Minimum Energy Path).
 * 
 * 	 | Finding saddle points on high dimensional surfaces.
 * 	 |    J. Chem. Phys., Vol. 111, No. 15, 15 October 1999.
 * 
 * 	@author msellers
 */

public class IntegratorDimerMEP extends IntegratorBox implements AgentSource {

	public ISimulation sim;
	public Box boxMin;
	public AtomAgentManager atomAgent0, atomAgentMin;
	public PotentialCalculationForceSum force0, forceMin;
	public IteratorDirective allatoms;
	public IRandom random;
	public MeterPotentialEnergy energyBox0, energyBoxMin;

	public IVectorRandom [] N, Nstar;
	public IVector NDelta, NstarDelta;
	public IVector workVector1;
	public IVector workVector2;
	public IVector [] saddle;
	public IVector [] F0, Fmin, Fmin2;
	public IVector [] THETA, THETAstar;
	public IVector [] Fperp, Fminperp, Fmin2perp;
	public IVector [] Fstar, Fminstar, Fmin2star, Fstarperp;
	public IVector [] Fpara;
	public Species [] movableSpecies;
	public AtomArrayList list, listMin;
	public int movableAtoms;
	public double energy;
	public double deltaR;
	public double deltaTheta;
	public double stepLength;
	public double dTheta;
	public double Frot, dFrot;
	public double Fprimerot;
	public double sinDtheta, cosDtheta;
	public int rotCounter, counter;
	public boolean rotate;
	
	
	public IntegratorDimerMEP(ISimulation sim, PotentialMaster potentialMaster, IVector [] asaddle, double deltaR, double stepLength, Species[] species) {
		this(sim, potentialMaster, sim.getRandom(), stepLength, deltaR, 1.0, asaddle, species);
	}
	
	public IntegratorDimerMEP(ISimulation aSim, PotentialMaster potentialMaster, IRandom arandom, double astepLength, double adeltaR, double temperature, IVector [] asaddle, Species[] aspecies) {
		super(potentialMaster, temperature);
		this.random = arandom;
		this.saddle = asaddle;
		this.sim = aSim;
		this.force0 = new PotentialCalculationForceSum();
		this.forceMin = new PotentialCalculationForceSum();
		this.allatoms = new IteratorDirective();
		this.stepLength = astepLength;
		this.deltaR = adeltaR;
		this.movableSpecies = aspecies;
		
		dTheta = 10E-4;
		dFrot = 0.1;
		rotCounter = 0;
		counter = 1;
		Frot = 1;
		rotate = true;
		
		sinDtheta = Math.sin(dTheta)*deltaR;
		cosDtheta = Math.cos(dTheta)*deltaR;
				
	}
	

	
	public void doStepInternal(){
		
		// Orient half-dimer on minimum energy path
		rotateDimerNewton();
		
		// Step half-dimer toward the local energy minimum
		walkDimer();
		
		System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2)
					+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(2));
		
		// Check and see if we're at the minimum energy
		if(counter%4==0){energyDimer();}
		
		counter++;
	}
	
	
	protected void setup() throws ConfigurationOverlapException{
		super.setup();
			       
	        movableAtoms = 0;
	        
	        for(int i=0; i<movableSpecies.length; i++){
	            movableAtoms += box.getMoleculeList(movableSpecies[i]).getAtomCount();
	        }
	        
	        N = new IVectorRandom [movableAtoms];
	        F0 = new IVector [movableAtoms];
	        Fmin = new IVector [movableAtoms];
	        Fmin2 = new IVector [movableAtoms];
	        Nstar = new IVectorRandom [movableAtoms];
	        THETA = new IVector [movableAtoms];
	        THETAstar = new IVector [movableAtoms];
	        Fperp = new IVector [movableAtoms];
	        Fminperp = new IVector [movableAtoms];
	        Fmin2perp = new IVector [movableAtoms];
	        Fstar = new IVector [movableAtoms];
	        Fminstar = new IVector [movableAtoms];
	        Fmin2star = new IVector [movableAtoms];
	        Fstarperp = new IVector [movableAtoms];
	        Fpara = new IVector [movableAtoms];
	        
	        for (int i=0; i<movableAtoms; i++){
	            
	            N[i] = (IVectorRandom)box.getSpace().makeVector();
	            Nstar[i] = (IVectorRandom)box.getSpace().makeVector();
	            F0[i] = box.getSpace().makeVector();
	            Fmin[i] = box.getSpace().makeVector();
	            Fmin2[i] = box.getSpace().makeVector();
	            Nstar[i] = (IVectorRandom)box.getSpace().makeVector();
	            NDelta = box.getSpace().makeVector();
	            NstarDelta = box.getSpace().makeVector();
	            THETA[i] = box.getSpace().makeVector();
	            THETAstar[i] = box.getSpace().makeVector();
	            Fperp[i] = box.getSpace().makeVector();
	            Fminperp[i] = box.getSpace().makeVector();
	            Fmin2perp[i] = box.getSpace().makeVector();
	            Fstar[i] = box.getSpace().makeVector();
	            Fminstar[i] = box.getSpace().makeVector();
	            Fmin2star[i] = box.getSpace().makeVector();
	            Fstarperp[i] = box.getSpace().makeVector();
	            Fpara[i] = box.getSpace().makeVector();
	        }
		
		// Use random unit array for N to generate dimer from saddle
		
		// Normalize N
		double mag = 0;
		for (int i=0; i<N.length; i++){
			N[i].setRandomSphere(random);
			mag += N[i].squared();
		}
		mag = 1.0 / Math.sqrt(mag);
		for (int i=0; i<N.length; i++){
			N[i].Ea1Tv1(mag, N[i]);
		}
		
		boxMin = new Box(box.getBoundary());
        sim.addBox(boxMin);
         
        energyBox0 = new MeterPotentialEnergy(this.potential);
        energyBoxMin = new MeterPotentialEnergy(this.potential);
         
        energyBox0.setBox(box);
        energyBoxMin.setBox(boxMin);
         
		// Offset Rmin (half-dimer end) from initial configuration, along N.		
		atomAgent0 = new AtomAgentManager(this, box);
		atomAgentMin = new AtomAgentManager(this, boxMin);
		
		force0.setAgentManager(atomAgent0);
		forceMin.setAgentManager(atomAgentMin);
		
		Species [] species = sim.getSpeciesManager().getSpecies();
		
		for(int i=0; i<species.length; i++){
			boxMin.setNMolecules(species[i], box.getNMolecules(species[i]));
		}
		
		// Setup atoms at saddle point, and create dimer for minimum path tracing
		for(int i=0; i<box.atomCount(); i++){
		//	((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition().E(saddle[i]);
			((IAtomPositioned)boxMin.getLeafList().getAtom(i)).getPosition().E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());	
		}
				
		// Atom list for movable and offset atoms
		list = new AtomArrayList();
        listMin = new AtomArrayList();
		
        for(int i=0; i<movableSpecies.length; i++){
            list.addAll(box.getMoleculeList(movableSpecies[i]));
            listMin.addAll(boxMin.getMoleculeList(movableSpecies[i]));
        }
		
        //Offset replica
        for(int i=0; i<N.length; i++){
            ((IAtomPositioned)listMin.getAtom(i)).getPosition().PEa1Tv1(deltaR, N[i]);
        }
		
        // Write out initial configuration
        System.out.println("----Dimer Minimum Energy Path");
        System.out.println(N[0].x(0)+"     "+N[0].x(1)+"     "+N[0].x(2));       
        System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2)
                +"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(2));
        
		// Calculate F's
		dimerForces(Fmin, F0, Fmin2);
		
	}
	
	public void rotateDimerNewton(){
        
	    if(Frot<dFrot){rotate=false;}
	    
		while(rotate){
			
			// Calculate F|_
			dimerForcePerp(N, Fmin, Fperp);
			
			// Find THETA and normalize
			double mag = 0;
			for(int i=0; i<Fperp.length; i++){
				THETA[i].E(Fperp[i]);
				mag += THETA[i].squared();
			}
			mag = 1.0 / Math.sqrt(mag);
			
			for(int i=0; i<THETA.length; i++){
				THETA[i].Ea1Tv1(mag, THETA[i]);
			}
			
			// Find R*'s
			for(int i=0; i<N.length; i++){
				IVector workVector;
				
				NstarDelta.Ea1Tv1(cosDtheta, N[i]);
				NstarDelta.PEa1Tv1(sinDtheta, THETA[i]);
							
				// Rmin*
				workVector = ((IAtomPositioned)listMin.getAtom(i)).getPosition();
				workVector.E(((IAtomPositioned)list.getAtom(i)).getPosition());
				workVector.PE(NstarDelta);
								
				// Find N*
				Nstar[i].Ea1Tv1(1.0/deltaR, NstarDelta);
			}
			
			// Calculate F*'s
			dimerForces(Fminstar, Fstar, Fmin2star);
			
			// Calculate F*|_
			dimerForcePerp(Nstar, Fminstar, Fstarperp);
			
			// Find THETA* and normalize
			mag = 0;
			for(int i=0; i<Fstarperp.length; i++){
				THETAstar[i].E(Fstarperp[i]);
				mag += THETAstar[i].squared();
			}
			mag = 1.0 / Math.sqrt(mag);
			
			for(int i=0; i<THETAstar.length; i++){
				THETAstar[i].Ea1Tv1(mag, THETAstar[i]);
			}
			
			// Compute scalar rotational force and change in rotational force over length dTheta
			Frot = 0;
			Fprimerot = 0;
			for(int i=0; i<Fperp.length; i++){
				Frot += Fperp[i].dot(THETA[i]);
				Frot += Fstarperp[i].dot(THETAstar[i]);
			}
			Frot = Frot/2.0;
			
			for(int i=0; i<Fstar.length; i++){
				Fprimerot += Fstarperp[i].dot(THETAstar[i]);
			}
			for(int i=0; i<Fstar.length; i++){
				Fprimerot -= Fperp[i].dot(THETA[i]);
			}
			
			Fprimerot = Fprimerot / dTheta;
			
			
			// Find actual rotation angle to minimize energy
			
			// Adatom on surface 
			deltaTheta = -0.5 * Math.atan2(2.0*Frot,Fprimerot) - dTheta/2.0;
			
			
			if(Fprimerot<0){deltaTheta = deltaTheta + Math.PI/2.0;}
			
			// Rotate dimer to new positions and set new N
			double sindeltaTheta = Math.sin(deltaTheta) * deltaR;
			double cosdeltaTheta = Math.cos(deltaTheta) * deltaR;
			
			for(int i=0; i<N.length; i++){
				IVector workVector;
				
				NDelta.Ea1Tv1(cosdeltaTheta, N[i]);
				NDelta.PEa1Tv1(sindeltaTheta, THETA[i]);
							
				// Rmin**
				workVector = ((IAtomPositioned)listMin.getAtom(i)).getPosition();
				workVector.E(((IAtomPositioned)list.getAtom(i)).getPosition());
				workVector.PE(NDelta);
			}
			
			// Calculate new F's
			dimerForces(Fmin, F0, Fmin2);
			
			// Calculate new Normal
			dimerNormal();
			
			rotCounter++;
			
			if(rotCounter>5){rotate=false;}
		}
		//Reset rotate
		rotate=true;
	}
	
	// Moves the half-dimer stepLength*N towards the energy minimum
	public void walkDimer(){
		
		IVector workvector;
		workvector = box.getSpace().makeVector();
		
		for(int i=0; i<N.length; i++){
			workvector.Ea1Tv1(stepLength, N[i]);
			((IAtomPositioned)list.getAtom(i)).getPosition().PE(workvector);
			((IAtomPositioned)listMin.getAtom(i)).getPosition().PE(workvector);
		}
		
		dimerForces(Fmin, F0, Fmin2);
		
		dimerNormal();
	}
	
	public void energyDimer(){
		
		double eMin, e0;
		e0 = energyBox0.getDataAsScalar();
		eMin = energyBoxMin.getDataAsScalar();
		
		System.out.println(e0);
		
		if(e0 < eMin){ 
			System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2)
					+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)listMin.getAtom(0)).getPosition().x(2));
						
			System.out.println("Box0 "+e0);
			System.out.println("BoxMin "+eMin);
			System.exit(1);
		}
		
	}
	
	// Compute Normal vector for dimer orientation
	protected void dimerNormal(){
	
		IVector workvector;
		workvector = box.getSpace().makeVector();
		
		// N =  (Rmin - R0) / (deltaR)
		for (int i=0; i<N.length; i++){	
		
			workvector.E(((IAtomPositioned)listMin.getAtom(i)).getPosition());
			workvector.ME(((IAtomPositioned)list.getAtom(i)).getPosition());
			
			N[i].Ea1Tv1(1.0/(deltaR), workvector);
		}
	}
		
	// Reset forces in boxes 0 and min, call calculate, and copy over new forces
	protected void dimerForces(IVector [] aF1, IVector [] aF, IVector [] aF2){
		force0.reset();
		forceMin.reset();
		
		potential.calculate(box, allatoms, force0);
		potential.calculate(boxMin, allatoms, forceMin);
		
		// Copy forces of dimer ends (R1, R2) to local array
		for(int i=0; i<aF1.length; i++){
			
			aF[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent0.getAgent(list.getAtom(i))).force());
			aF1[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgentMin.getAgent(listMin.getAtom(i))).force());
			aF2[i].Ea1Tv1(2.0, aF[i]);
			aF2[i].ME(aF1[i]);
			
		}
	}
	
	// Calculate the force, aF, perpendicular to dimer orientation N.
	protected void dimerForcePerp(IVectorRandom [] aN, IVector [] aF1, IVector [] aFperp){
		
		double mag1 = 0;
		
		// F|_ = F1 - (F1[dot]N)*N
		for(int i=0; i<aF1.length; i++){
			mag1 += aF1[i].dot(aN[i]);
		}
		for (int i=0; i<aFperp.length; i++){
			aFperp[i].E(aF1[i]);
			aFperp[i].PEa1Tv1(-mag1, aN[i]);	
			aFperp[i].TE(1.0/deltaR);
		}

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
