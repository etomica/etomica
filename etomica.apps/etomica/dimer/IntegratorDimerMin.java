package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.api.IAtom;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IRandom;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IVector;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.IVectorRandom;
import etomica.space.Space;
import etomica.units.ElectronVolt;


/**
 * 	Henkelman's Dimer Method (Minimum Energy Path).
 * 
 * 	 | Finding saddle points on high dimensional surfaces.
 * 	 |    J. Chem. Phys., Vol. 111, No. 15, 15 October 1999.
 * 
 * 	@author msellers
 */

public class IntegratorDimerMin extends IntegratorBox implements AgentSource {

	public ISimulation sim;
	public IBox boxMin;
	public AtomAgentManager atomAgent0, atomAgentMin;
	public PotentialCalculationForceSum force0, forceMin;
	public IteratorDirective allatoms;
	public IRandom random;
	public FileWriter fileWriter;
	public MeterPotentialEnergy energyBox0, energyBoxMin;
	public ActivityIntegrate activityIntegrate;

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
	public ISpecies [] movableSpecies;
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
	public double e0prev;
	public int rotCounter, counter;
	public boolean rotate, normalD;
	public String file;
	public WriteConfiguration writer;
	private final Space space;
	
	
	public IntegratorDimerMin(ISimulation sim, PotentialMaster potentialMaster,
			                  ISpecies[] species, String fileName,
			                  Boolean normalDir, Space _space) {
		this(sim, potentialMaster, sim.getRandom(), 1.0, species, fileName, normalDir, _space);
	}
	
	public IntegratorDimerMin(ISimulation aSim, PotentialMaster potentialMaster,
			                  IRandom arandom, double temperature,
			                  ISpecies[] aspecies, String fileName,
			                  Boolean normalDir, Space _space) {
		super(potentialMaster, temperature);
		this.random = arandom;
		this.sim = aSim;
		this.force0 = new PotentialCalculationForceSum();
		this.forceMin = new PotentialCalculationForceSum();
		this.allatoms = new IteratorDirective();
		this.movableSpecies = aspecies;
		this.file = fileName;
		this.normalD = normalDir;
		this.space = _space;
		
		stepLength = 0.05;
		deltaR = 10E-4;
		dTheta = 10E-4;
		dFrot = 0.1;
		rotCounter = 0;
		counter = 1;
		Frot = 1;
		rotate = true;		
		e0prev = 0;
	}
		
	public void doStepInternal(){
		
		// Orient half-dimer on minimum energy path
		rotateDimerNewton();
		
		// Write energy to file
        try{
            fileWriter = new FileWriter(file+"_minimum_path");
            fileWriter.write(ElectronVolt.UNIT.fromSim(energyBox0.getDataAsScalar())+"\n");
        }catch(IOException e) {
          
        }
		
        // Check and see if we're at the minimum energy
        if(counter>100){
        energyDimer();
        }
        
		// Step half-dimer toward the local energy minimum
		walkDimer();
		
		System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2));       
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
	            N[i] = (IVectorRandom)space.makeVector();
	            Nstar[i] = (IVectorRandom)space.makeVector();
	            F0[i] = space.makeVector();
	            Fmin[i] = space.makeVector();
	            Fmin2[i] = space.makeVector();
	            Nstar[i] = (IVectorRandom)space.makeVector();
	            NDelta = space.makeVector();
	            NstarDelta = space.makeVector();
	            THETA[i] = space.makeVector();
	            THETAstar[i] = space.makeVector();
	            Fperp[i] = space.makeVector();
	            Fminperp[i] = space.makeVector();
	            Fmin2perp[i] = space.makeVector();
	            Fstar[i] = space.makeVector();
	            Fminstar[i] = space.makeVector();
	            Fmin2star[i] = space.makeVector();
	            Fstarperp[i] = space.makeVector();
	            Fpara[i] = space.makeVector();
	        }  
	        
		boxMin = new Box(box.getBoundary(), space);
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
		
		ISpecies [] species = sim.getSpeciesManager().getSpecies();
		
		for(int i=0; i<species.length; i++){
			boxMin.setNMolecules(species[i], box.getNMolecules(species[i]));
		}
		
		// Read in coordinates for boxMin atom locations
		ConfigurationFile configFile = new ConfigurationFile(file+"_fine_A_saddle");
    	configFile.initializeCoordinates(boxMin);
    	writer = new WriteConfiguration(space);
    	writer.setConfName(file+"_A_minimum");
    	
    	if(normalD==true){
    		// Read in coordinates for opposite boxMin atom locations
    		ConfigurationFile configFile1 = new ConfigurationFile(file+"_fine_B_saddle");
        	configFile1.initializeCoordinates(boxMin);
        	writer.setConfName(file+"_B_minimum");
    	}
    		
		// Atom list for movable and offset atoms
		list = new AtomArrayList();
        listMin = new AtomArrayList();
		
		for(int i=0; i<movableSpecies.length; i++){
            IAtomSet molecules = box.getMoleculeList(movableSpecies[i]);
            IAtomSet molecules1 = boxMin.getMoleculeList(movableSpecies[i]);
            for (int j=0; j<molecules.getAtomCount(); j++) {
                list.add(((IMolecule)molecules.getAtom(j)).getChildList().getAtom(0));
                listMin.add(((IMolecule)molecules1.getAtom(j)).getChildList().getAtom(0));
            }
		}  
		
        dimerNormal();
        
        // Write out initial configuration
        System.out.println(file+" ***Dimer Minima Search***");
        System.out.println(N[0].x(0)+"     "+N[0].x(1)+"     "+N[0].x(2));       
        System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2));
        
		// Calculate F's
		dimerForces(Fmin, F0, Fmin2);
		
	}
	
	public void rotateDimerNewton(){
        
	    dTheta = 10E-5;
	    deltaTheta = 1.0;
	    Frot = 1.0;
	    
		while(true){
			
			// Calculate F|_
			dimerForcePerp(N, Fmin, Fperp);
			
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
			
            // NEWTON'S LINE MINIMIZATION FOR ROTATIONAL FORCE	
            Frot = 0;
            for(int i=0; i<Fperp.length; i++){
                  Frot += Fperp[i].dot(THETA[i]);
            } 
            // Leave loop if Frot at current location is small
            if(Frot<dFrot){break;}
                        
            sinDtheta = Math.sin(dTheta);
            cosDtheta = Math.cos(dTheta);
            
			// Find Nstar and THETAstar after dTheta rotation
	        IVector workVectorN1 = space.makeVector();
			for(int i=0; i<N.length; i++){
				workVectorN1.Ea1Tv1(cosDtheta, N[i]);
				workVectorN1.PEa1Tv1(sinDtheta, THETA[i]);
				Nstar[i].E(workVectorN1);
				workVectorN1.Ea1Tv1(-sinDtheta, N[i]);
				workVectorN1.PEa1Tv1(cosDtheta, THETA[i]);
				THETAstar[i].E(workVectorN1);
			}
                        
            // Use Nstar to offset(rotate) replicas
            IVector workVector1 = space.makeVector();
            for(int i=0; i<Nstar.length; i++){
                workVector1.E(((IAtomPositioned)list.getAtom(i)).getPosition());
                workVector1.PEa1Tv1(deltaR, Nstar[i]);
                ((IAtomPositioned)listMin.getAtom(i)).getPosition().E(workVector1);
            }
			
			// Calculate F*'s
			dimerForcesStar(Fminstar, Fmin2star, F0);
			
			// Calculate F*|_
			dimerForcePerp(Nstar, Fminstar, Fstarperp);
			
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
			deltaTheta = -0.5 * Math.atan(2.0*Frot/Fprimerot) - dTheta/2.0;
			if(Fprimerot>0){deltaTheta = deltaTheta + Math.PI/2.0;}
			System.out.println("Frot "+Frot+"    Fprimerot "+Fprimerot);
			
             // Check deltaTheta vs. dTheta and adjust step size
            if (deltaTheta < 0){
                    dTheta /= 10;
                    System.out.println("dTheta changed to "+dTheta);
                    if (deltaTheta < -dTheta){
                        deltaTheta = -dTheta;
                    }
            }
            
            double sindeltaTheta = Math.sin(deltaTheta);
            double cosdeltaTheta = Math.cos(deltaTheta);
            
            // Find N**
            IVector workVectorN2 = space.makeVector();
            for(int i=0; i<N.length; i++){               
                workVectorN2.Ea1Tv1(cosdeltaTheta, Nstar[i]);
                workVectorN2.PEa1Tv1(sindeltaTheta, THETAstar[i]);
                N[i].E(workVectorN2);
            }
            
            // Use new N to offset(rotate) replica
            IVector workVector2 = space.makeVector();
            for(int i=0; i<N.length; i++){             
                workVector2.E(((IAtomPositioned)list.getAtom(i)).getPosition());
                workVector2.PEa1Tv1(deltaR, N[i]);
                ((IAtomPositioned)listMin.getAtom(i)).getPosition().E(workVector2);
            }     
			
			// Calculate new F's
			dimerForces(Fmin, F0, Fmin2);
			
			// Calculate new Normal
			dimerNormal();
			
			rotCounter++;
			
			//if(rotCounter>15){
	        //        break;
			//}
		}
	}
	
	// Moves the half-dimer stepLength*N towards the energy minimum
	public void walkDimer(){
		
		IVector workvector;
		workvector = space.makeVector();
		
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
		
		
    	try { 
            fileWriter.write(ElectronVolt.UNIT.fromSim(e0)+"\n");
            }catch(IOException e) {
                System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            }
		
		
		System.out.println(ElectronVolt.UNIT.fromSim(e0));
		
		if(e0>e0prev){ 
		    System.out.println(file+" +++Dimer Minimum Found+++");
			System.out.println("Box0 = "+ElectronVolt.UNIT.fromSim(e0)+"eV");
			eMin = energyBoxMin.getDataAsScalar();
			System.out.println("BoxMin = "+ElectronVolt.UNIT.fromSim(eMin)+"eV");
			
			try { 
	            fileWriter.close();
	            }catch(IOException e) {
	                System.err.println("Cannot open file, caught IOException: " + e.getMessage());
	            }
	            
            writer.setBox(box);
            writer.actionPerformed();
            activityIntegrate.setMaxSteps(0);
		}
		e0prev = e0;
	}
	
	// Compute Normal vector for dimer orientation
	protected void dimerNormal(){
        double mag=0;
        IVector workvector;
        workvector = space.makeVector();
        
        // N =  (Rmin - R0) / (deltaR)
        for (int i=0; i<N.length; i++){ 
            workvector.E(((IAtomPositioned)listMin.getAtom(i)).getPosition());
            workvector.ME(((IAtomPositioned)list.getAtom(i)).getPosition());
            N[i].E(workvector);
            mag += workvector.squared();
        }
        
        mag = 1.0/Math.sqrt(mag);
        
        for (int i=0; i<N.length; i++){
            N[i].TE(mag);
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
	
	protected void dimerForcesStar(IVector [] aF1star, IVector [] aF2star, IVector [] aF){
        forceMin.reset();
        potential.calculate(boxMin, allatoms, forceMin);
        
     // Copy forces of dimer end and center (R1, R) to local array
        for(int i=0; i<aF1star.length; i++){
            aF1star[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgentMin.getAgent(listMin.getAtom(i))).force());
            aF2star[i].Ea1Tv1(2.0, aF[i]);
            aF2star[i].ME(aF1star[i]);
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
	
	public void setActivityIntegrate(ActivityIntegrate ai){
		activityIntegrate = ai;
	}
	
	public Class getAgentClass() {
		return IntegratorVelocityVerlet.MyAgent.class;
	}

	public Object makeAgent(IAtom a) {
		return new IntegratorVelocityVerlet.MyAgent(space);
	}

	public void releaseAgent(Object agent, IAtom atom) {
		// TODO Auto-generated method stub	
	}
	
}
