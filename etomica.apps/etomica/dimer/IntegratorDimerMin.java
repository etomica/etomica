package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
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
import etomica.simulation.ISimulation;
import etomica.space.IVector;
import etomica.space.IVectorRandom;
import etomica.species.Species;
import etomica.units.ElectronVolt;
import etomica.util.IRandom;

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
	public Box boxMin;
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
	public boolean rotate, normalD;
	public String file;
	public WriteConfiguration writer;
	
	
	public IntegratorDimerMin(ISimulation sim, PotentialMaster potentialMaster, Species[] species, String fileName, Boolean normalDir) {
		this(sim, potentialMaster, sim.getRandom(), 1.0, species, fileName, normalDir);
	}
	
	public IntegratorDimerMin(ISimulation aSim, PotentialMaster potentialMaster, IRandom arandom, double temperature, Species[] aspecies, String fileName, Boolean normalDir) {
		super(potentialMaster, temperature);
		this.random = arandom;
		this.sim = aSim;
		this.force0 = new PotentialCalculationForceSum();
		this.forceMin = new PotentialCalculationForceSum();
		this.allatoms = new IteratorDirective();
		this.movableSpecies = aspecies;
		this.file = fileName;
		this.normalD = normalDir;
		
		stepLength = 10E-3;
		deltaR = 10E-4;
		dTheta = 10E-4;
		dFrot = 0.01;
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
		
		System.out.println(((IAtomPositioned)list.getAtom(0)).getPosition().x(0)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(1)+"     "+((IAtomPositioned)list.getAtom(0)).getPosition().x(2));
		
		// Write energy to file
        try{
            fileWriter = new FileWriter(file+"_minimum_path");
            fileWriter.write(ElectronVolt.UNIT.fromSim(energyBox0.getDataAsScalar())+"\n");
            fileWriter.close();
        }catch(IOException e) {
          
        }
		
		// Check and see if we're at the minimum energy
		energyDimer();
		
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
		
		// Read in coordinates for boxMin atom locations
		ConfigurationFile configFile = new ConfigurationFile(file+"_1_saddle");
    	configFile.initializeCoordinates(boxMin);
    	writer = new WriteConfiguration();
    	writer.setConfName(file+"_A_minimum");
    	
    	if(normalD==true){
    		// Read in coordinates for opposite boxMin atom locations
    		ConfigurationFile configFile1 = new ConfigurationFile(file+"_2_saddle");
        	configFile1.initializeCoordinates(boxMin);
        	writer.setConfName(file+"_B_minimum");
    	}
    		
		// Atom list for movable and offset atoms
		list = new AtomArrayList();
        listMin = new AtomArrayList();
		
        for(int i=0; i<movableSpecies.length; i++){
            list.addAll(box.getMoleculeList(movableSpecies[i]));
            listMin.addAll(boxMin.getMoleculeList(movableSpecies[i]));
        }     
		
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
                ((IAtomPositioned)listMin.getAtom(i)).getPosition().E(workVector1);
            }
			
			// Calculate F*'s
			dimerForcesStar(F0, Fminstar, Fstar, Fmin2star);
			
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
			
			// Adatom on surface 
			deltaTheta = -0.5 * Math.atan(2.0*Frot/Fprimerot) - dTheta/2.0;
			
			
			if(Fprimerot>0){deltaTheta = deltaTheta + Math.PI/2.0;}
			
            //System.out.println(energy+"    "+Frot+"     "+Fprimerot+"     "+deltaTheta);
            
             // Check deltaTheta vs. dTheta and adjust step size
            if (deltaTheta < 0){
                    dTheta /= 100;
                    System.out.println("dTheta changed to "+dTheta);
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
            
            // Use new N to offset(rotate) replica
            IVector workVector2 = box.getSpace().makeVector();
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
		
    	try { 
            fileWriter.write(ElectronVolt.UNIT.fromSim(e0)+"\n");
            }catch(IOException e) {
                System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            }
		
		
		System.out.println(ElectronVolt.UNIT.fromSim(e0));
		
		if(e0 < eMin){ 
		    System.out.println(file+" +++Dimer Minimum Found+++");
			System.out.println("Box0 = "+ElectronVolt.UNIT.fromSim(e0));
			System.out.println("BoxMin = "+ElectronVolt.UNIT.fromSim(eMin));
			
			try { 
	            fileWriter.close();
	            }catch(IOException e) {
	                System.err.println("Cannot open file, caught IOException: " + e.getMessage());
	            }
	            
            writer.setBox(box);
            writer.actionPerformed();
            activityIntegrate.setMaxSteps(0);
		}
		
	}
	
	// Compute Normal vector for dimer orientation
	protected void dimerNormal(){
        double mag=0;
        IVector workvector;
        workvector = box.getSpace().makeVector();
        
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
	
	protected void dimerForcesStar(IVector [] aF, IVector [] aF1star, IVector [] aF2star, IVector []aFstar){
        forceMin.reset();
        potential.calculate(boxMin, allatoms, forceMin);
        
     // Copy forces of dimer end and center (R1, R) to local array
        for(int i=0; i<aF1star.length; i++){
            aF1star[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgentMin.getAgent(listMin.getAtom(i))).force());
            aFstar[i].E(aF[i]);
            aF2star[i].Ea1Tv1(2.0, aFstar[i]);
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
		return new IntegratorVelocityVerlet.MyAgent(box.getSpace());
	}

	public void releaseAgent(Object agent, IAtom atom) {
		// TODO Auto-generated method stub	
	}
	
}
