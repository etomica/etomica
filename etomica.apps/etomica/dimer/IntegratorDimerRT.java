package etomica.dimer;

import etomica.atom.AtomAgentManager;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVector;
import etomica.space3d.Vector3D;
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
	public IVector [] N, Nstar;
	public IVector NDelta, NstarDelta;
	public double deltaR;
	public double dTheta;
	public double deltaTheta;
	public IVector [] R1star, R2star, THETA, THETAstar, THETAstar2;
	public IVector [] F, F1, F2;
	public IVector [] Fperp, F1perp, F2perp;
	public IVector [] Fstar, F1star, F2star, Fstarperp;
	public double Frot, Fprimerot;
	public IRandom random1;
	public PotentialCalculationForceSum force1, force2;
	public AtomAgentManager atomAgent1, atomAgent2;
	public IteratorDirective allatoms;
	public double sinDtheta, cosDtheta;
	
	public IntegratorDimerRT(ISimulation sim, PotentialMaster potentialMaster) {
		this(sim, potentialMaster, sim.getRandom(), 0.05, 1.0);
	}
	
	public IntegratorDimerRT(ISimulation aSim, PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature) {
		super(potentialMaster, temperature);
		this.random1 = random;
		this.sim = aSim;
		this.force1 = new PotentialCalculationForceSum();
		this.force2 = new PotentialCalculationForceSum();
		this.allatoms = new IteratorDirective();
		
		deltaR = 1E-2;
		dTheta = 1E-4;
		sinDtheta = Math.sin(dTheta)*deltaR;
		cosDtheta = Math.cos(dTheta)*deltaR;
	}
	
	
	public void doStepInternal(){
		
		
	
		
	}
	
		
	// Takes in a current configuration of atoms (Rc) and creates a dimer of their positions (R1 and R2).
	// (R2)-------[Rc]-------(R1) ===> N

	protected void createDimer(){
			
		//Use random unit array for N to generate dimer.
		N = new IVector [box.atomCount()];
		Vector3D n = new Vector3D();
		
		// Normalize N
		double mag = 0;
		for (int i=0; i<N.length; i++){ 
			n.setRandomSphere(random1);
			N[i].E(n);
			mag += N[i].squared();
		}
		for (int i=0; i<N.length; i++){
			N[i].Ea1Tv1(1.0/Math.sqrt(mag), N[i]);
		}
		mag = 0;
		
		//Offset R1 and R2 (dimer ends) from initial configuration, along N.
		box1 = new Box(box.getBoundary());
		box2 = new Box(box.getBoundary());
		
		sim.addBox(box1);
		sim.addBox(box2);
		
		atomAgent1 = new AtomAgentManager(this, box1);
		atomAgent2 = new AtomAgentManager(this, box2);
		
		force1.setAgentManager(atomAgent1);
		force2.setAgentManager(atomAgent2);
		
		AtomSet list = box.getLeafList();
		
		Species [] species = sim.getSpeciesManager().getSpecies();
		
		for(int i=0; i<species.length; i++){
			box1.setNMolecules(species[i], box.getNMolecules(species[i]));
			box2.setNMolecules(species[i], box.getNMolecules(species[i]));
		}
		
		AtomSet list1 = box1.getLeafList();
		AtomSet list2 = box2.getLeafList();

		for(int i=0; i<box.atomCount(); i++){
			((IAtomPositioned)list1.getAtom(i)).getPosition().E(((IAtomPositioned)list.getAtom(i)).getPosition());
			((IAtomPositioned)list1.getAtom(i)).getPosition().PEa1Tv1(deltaR, N[i]);
			
			((IAtomPositioned)list2.getAtom(i)).getPosition().E(((IAtomPositioned)list.getAtom(i)).getPosition());
			((IAtomPositioned)list2.getAtom(i)).getPosition().PEa1Tv1(-deltaR, N[i]);
		}
		
	}
	
	// Rotate the dimer to align it with the lowest curvature mode of the potential energy surface, Newton style.
	public void rotateDimerNewton(){
		
		dimerForce();
		
		//Copy out forces of dimer ends to local array
		for(int i=0; i<box.atomCount(); i++){
			F1[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent1.getAgent(box1.getLeafList().getAtom(i))).force());
			F2[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent1.getAgent(box1.getLeafList().getAtom(i))).force());
			F[i].Ev1Mv2(F1[i], F2[i]);
		}
		
		dimerForcePerp(N, F1, F2, Fperp);
		
		double mag = 0;
		for(int i=0; i<Fperp.length; i++){
			THETA[i].E(Fperp[i]);
			mag += THETA[i].squared();
		}
		for(int i=0; i<THETA.length; i++){
			THETA[i].Ea1Tv1(1.0/Math.sqrt(mag), THETA[i]);
		}
				
		for(int i=0; i<box.atomCount(); i++){
			IVector workVector;
			
			NstarDelta.Ea1Tv1(cosDtheta, N[i]);
			NstarDelta.PEa1Tv1(sinDtheta, THETA[i]);
						
			// R1*
			workVector = ((IAtomPositioned)box1.getLeafList().getAtom(i)).getPosition();
			workVector.E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
			workVector.PE(NstarDelta);
			
			// R2*
			workVector = ((IAtomPositioned)box2.getLeafList().getAtom(i)).getPosition();
			workVector.E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
			workVector.ME(NstarDelta);
		}
		
		// Calculate F*'s
		dimerForce();
		
		//Copy out forces of dimer ends to local array
		for(int i=0; i<box.atomCount(); i++){
			F1star[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent1.getAgent(box1.getLeafList().getAtom(i))).force());
			F2star[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent1.getAgent(box1.getLeafList().getAtom(i))).force());
			Fstar[i].Ev1Mv2(F1[i], F2[i]);
		}
		
		//Find Nstar
		mag = 0;
		for(int i=0; i<Fstar.length; i++){
			Nstar[i].E(Fstar[i]);
			mag += Nstar[i].squared();
		}
		for(int i=0; i<Nstar.length; i++){
			Nstar[i].Ea1Tv1(1/Math.sqrt(mag), Nstar[i]);
		}
		
		dimerForcePerp(Nstar, F1star, F2star, Fstarperp);
		
		//Find Theta*
		mag = 0;
		for(int i=0; i<Fstarperp.length; i++){
			THETAstar[i].E(Fstarperp[i]);
			mag += THETAstar[i].squared();
		}
		for(int i=0; i<THETAstar.length; i++){
			THETAstar[i].Ea1Tv1(1.0/Math.sqrt(mag), THETAstar[i]);
		}
		
		//Compute scalar rotational force and change in force
		Frot = 0;
		Fprimerot = 0;
		for(int i=0; i<Fperp.length; i++){
			Frot += Fperp[i].dot(THETA[i]);
		}
		for(int i=0; i<Fstar.length; i++){
			Fprimerot += Fstar[i].dot(THETAstar[i]);
			Fprimerot -= F[i].dot(THETA[i]);
		}
		Fprimerot = Fprimerot / dTheta;
		
		//Find actual rotation angle to minimize energy
		deltaTheta = -0.5 * Math.atan(2.0*Frot/Fprimerot);
		
		//Rotate dimer to new positions and set new N
		double sindeltaTheta = Math.sin(deltaTheta) * deltaR;
		double cosdeltaTheta = Math.cos(deltaTheta) * deltaR;
		
		for(int i=0; i<box.atomCount(); i++){
			IVector workVector;
			
			NDelta.Ea1Tv1(cosdeltaTheta, N[i]);
			NDelta.PEa1Tv1(sindeltaTheta, THETA[i]);
						
			// R1**
			workVector = ((IAtomPositioned)box1.getLeafList().getAtom(i)).getPosition();
			workVector.E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
			workVector.PE(NDelta);
			
			// R2**
			workVector = ((IAtomPositioned)box2.getLeafList().getAtom(i)).getPosition();
			workVector.E(((IAtomPositioned)box.getLeafList().getAtom(i)).getPosition());
			workVector.ME(NDelta);
			
			N[i].Ea1Tv1(1.0/deltaR, NDelta);
		}
	}
	
	
	
	
	protected void dimerForce(){
		force1.reset();
		force2.reset();
		
		potential.calculate(box1, allatoms, force1);
		potential.calculate(box2, allatoms, force2);
	}
	
	protected static void dimerForcePerp(IVector [] aN, IVector [] aF1, IVector [] aF2, IVector [] aF){
		
		double mag1 = 0;
		double mag2 = 0;
		
		// F = F1 - F1[dot]N * N  - F2  +  F2[dot]N * N
		for(int i=0; i<aF1.length; i++){
			mag1 += aF1[i].dot(aN[i]);
			mag2 += aF2[i].dot(aN[i]);
		}
		for (int i=0; i<aF.length; i++){
			aF[i].E(aF1[i]);
			aF[i].Ea1Tv1(-mag1, aN[i]);	
			aF[i].ME(aF2[i]);
			aF[i].PEa1Tv1(mag2, aN[i]);
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
