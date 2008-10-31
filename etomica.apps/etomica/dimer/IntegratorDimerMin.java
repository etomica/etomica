package etomica.dimer;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.CalcVibrationalModes;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorBox;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.PotentialCalculationForceSum;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;
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
	public AtomLeafAgentManager atomAgent0, atomAgentMin;
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
	public double eMin, e0;
	public int rotCounter, counter;
	public boolean rotate, normalD, minFound;
	public String file;
	public WriteConfiguration writer;
	private final ISpace space;
	public CalcVibrationalModes vib;
	
	
	public IntegratorDimerMin(ISimulation sim, IPotentialMaster potentialMaster,
			                  ISpecies[] species,
			                  Boolean normalDir, ISpace _space) {
		this(sim, potentialMaster, sim.getRandom(), 1.0, species, normalDir, _space);
	}
	
	public IntegratorDimerMin(ISimulation aSim, IPotentialMaster potentialMaster,
			                  IRandom arandom, double temperature,
			                  ISpecies[] aspecies, Boolean normalDir, ISpace _space) {
		super(potentialMaster, temperature);
		this.random = arandom;
		this.sim = aSim;
		this.force0 = new PotentialCalculationForceSum();
		this.forceMin = new PotentialCalculationForceSum();
		this.allatoms = new IteratorDirective();
		this.movableSpecies = aspecies;
		this.normalD = normalDir;
		this.space = _space;
		
		stepLength = 0.001;
		deltaR = 1E-3;
		dTheta = 1E-4;
		dFrot = 0.1;
		rotCounter = 0;
		counter = 0;
		Frot = 1;
		rotate = true;
		minFound = false;
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
		
	public void reset(){
	        stepLength = 0.005;
	        deltaR = 1E-3;
	        dTheta = 1E-4;
	        dFrot = 0.1;
	        rotCounter = 0;
	        counter = 0;
	        Frot = 1;
	        rotate = true;
	        minFound = false;
	    
	}
	
	public void initializeDimer(){
        ConfigurationFile config = new ConfigurationFile(file+"_saddle");
        config.initializeCoordinates(box);        
        
        // Read in coordinates for boxMin atom locations
        ConfigurationFile configFile = new ConfigurationFile(file+"_A_saddle");
        configFile.initializeCoordinates(boxMin);
        writer = new WriteConfiguration(space);
        writer.setConfName(file+"_A_minimum");
                
        if(normalD==true){
            // Read in coordinates for opposite boxMin atom locations
            ConfigurationFile configFile1 = new ConfigurationFile(file+"_B_saddle");
            configFile1.initializeCoordinates(boxMin);
            writer.setConfName(file+"_B_minimum");
        }
        
        try{
            fileWriter = new FileWriter(writer.getConfName()+"_path");
        }catch(IOException e) {
            
        }
        
        dimerNormal();
	}
	public void doStepInternal(){
	    // Step half-dimer toward the local energy minimum
	    walkDimer();
	    // Orient half-dimer on minimum energy path
        rotateDimerNewton();
	    	    
	    e0 = ElectronVolt.UNIT.fromSim(energyBox0.getDataAsScalar());
	    eMin = ElectronVolt.UNIT.fromSim(energyBoxMin.getDataAsScalar());
	    
        // Write energy to file
        try{
            fileWriter.write(e0+"\n");
        }catch(IOException e) {
       
        }
        
        /*
        dimerForces(Fmin, F0, Fmin2);
        
        if(counter>10){
            //Check slope of energy after step
            double slope=0;
            for(int i=0; i<F0.length; i++){
                slope += F0[i].dot(N[i]);
            }
            if(slope<0){
                quitSearch();
            }
        }
        
        //Find unit vector of avg. force
        double magAvg = 0;
        for(int i=0; i<N.length; i++){
            N[i].Ev1Pv2(F0[i], Fmin[i]);
            N[i].TE(0.5);
            magAvg+=N[i].squared();
        }
        magAvg=1.0/Math.sqrt(magAvg);
        for(int i=0; i<N.length; i++){
            N[i].TE(magAvg);
        }
        */
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
        
        if(potential instanceof PotentialMasterListDimer){
            this.addNonintervalListener(((PotentialMasterList)potential).getNeighborManager(boxMin));
            this.addIntervalAction(((PotentialMasterList)potential).getNeighborManager(boxMin)); 
         }
        
        energyBox0 = new MeterPotentialEnergy(this.potential);
        energyBoxMin = new MeterPotentialEnergy(this.potential);
         
        energyBox0.setBox(box);
        energyBoxMin.setBox(boxMin);
         
		// Offset Rmin (half-dimer end) from initial configuration, along N.		
		atomAgent0 = new AtomLeafAgentManager(this, box);
		atomAgentMin = new AtomLeafAgentManager(this, boxMin);
		
		force0.setAgentManager(atomAgent0);
		forceMin.setAgentManager(atomAgentMin);
		
//		ISpecies [] species = sim.getSpeciesManager().getSpecies();
		
		for(int i=0; i<sim.getSpeciesManager().getSpeciesCount(); i++){
			boxMin.setNMolecules(sim.getSpeciesManager().getSpecies(i),
					box.getNMolecules(sim.getSpeciesManager().getSpecies(i)));
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
	}
	
	/**
	 * Rotates the dimer in potential energy hyperspace and aligns it with the surface's lowest 
	 * curvature mode.  Using Newton's finite difference method, the rotational force perpendicular to the dimer's
	 * orientation is minimized.  Minimization criteria is dependent on dFrot and rotCounter.
	 * For minimum energy path tracing, rotCounter criteria is commented out to ensure convergence
	 * to lowest curvature mode before stepping the dimer.
	 */
	public void rotateDimerNewton(){
	    dimerNormal();
	    dimerForces(Fmin, F0, Fmin2);

	    Frot = 1.0;
	    
	    if(counter>10){
            //Check slope of energy after step
            double slope=0;
            for(int i=0; i<F0.length; i++){
                slope += F0[i].dot(N[i]);
            }
            if(slope<0){	
            	if(eMin>e0){
            		quitSearch();
            	}
           }
	    }
	    
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
			//System.out.println("Frot "+Frot+"    Fprimerot "+Fprimerot);
			
			/*
             // Check deltaTheta vs. dTheta and adjust step size
            if (deltaTheta < 0){
                    dTheta /= 10;
                    System.out.println("dTheta changed to "+dTheta);
                    if (deltaTheta < -dTheta){
                        deltaTheta = -dTheta;
                    }
            }
            */
			
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
						
			rotCounter++;
			
			if(rotCounter>100){
	                break;
			}
		}
	}
	
	/**
	 *  Moves the half-dimer stepLength*N towards the energy minimum.
	 */
	public void walkDimer(){
	    //System.out.println(e0);
		IVector workvector;
		workvector = space.makeVector();
		
		for(int i=0; i<N.length; i++){
			workvector.Ea1Tv1(stepLength, N[i]);
			((IAtomPositioned)list.getAtom(i)).getPosition().PE(workvector);
			((IAtomPositioned)listMin.getAtom(i)).getPosition().PE(workvector);
		}
		
	}
	
	/**
	 * Called after criteria is met in the rotateDimer method.  If F0dotN is 
	 * negative, the dimer has reached the minimum in an energy basin.
	 * 
	 */
	public void quitSearch(){
	    /*
	    double eMin;
	    eMin = ElectronVolt.UNIT.fromSim(energyBoxMin.getDataAsScalar());
	    System.out.println(file+" +++Dimer Minimum Found+++");
		System.out.println("Box0   = "+e0+" eV");
		System.out.println("BoxMin = "+eMin+" eV");
		*/
	    
	    minFound = true;
	    
        vib = new CalcVibrationalModes();
        vib.setup(box, super.potential, (IAtomSet)box.getMoleculeList(movableSpecies[0]), space);
        vib.actionPerformed();
        
		try { 
            fileWriter.close();
            }catch(IOException e) {
                System.err.println("Cannot open file, caught IOException: " + e.getMessage());
            }
        writer.setBox(box);
        writer.actionPerformed();
        //activityIntegrate.setMaxSteps(0);
	}
	
	/**
	 * Computes an n-dimensional vector in potential energy hyperspace, pointing 
	 * in the direction of the dimer (intersecting both simulation boxes).
	 * 
	 */
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
			
	/**
	 * Resets forces in boxes 0 and min, calls potential.calculate, and copies over
	 * the new forces to local arrays.  Approximates the force in virtual box 2.
	 * 
	 * @param aF1 Array of vectors holding the forces in box min
	 * @param aF Array of vectors holding the forces in box 0
	 * @param aF2 Array of vectors holding the forces in virtual box 2
	 */
	protected void dimerForces(IVector [] aF1, IVector [] aF, IVector [] aF2){
		force0.reset();
		forceMin.reset();
		
		potential.calculate(box, allatoms, force0);
		potential.calculate(boxMin, allatoms, forceMin);
		
		// Copy forces of dimer ends (R1, R2) to local array
		for(int i=0; i<aF1.length; i++){
			
			aF[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgent0.getAgent((IAtomLeaf)list.getAtom(i))).force());
			aF1[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgentMin.getAgent((IAtomLeaf)listMin.getAtom(i))).force());
			aF2[i].Ea1Tv1(2.0, aF[i]);
			aF2[i].ME(aF1[i]);
			
		}
	}
		
	/**
	 * Used in the middle of the rotateDimer method.  Resets forces in box min only (we are rotating
	 * the dimer around a fixed center).  Calls potential.calculate, and copies over the new forces
	 * to local array.  Approximates the force in virtual box 2.
	 * 
	 * @param aF1star Array of vectors holding the forces in box min
	 * @param aF2star Array of vectors holding the forces in virtual box 2
	 * @param aF Array of vectors holding the forces in box 0
	 */
	protected void dimerForcesStar(IVector [] aF1star, IVector [] aF2star, IVector [] aF){
        forceMin.reset();
        potential.calculate(boxMin, allatoms, forceMin);
        
     // Copy forces of dimer end and center (R1, R) to local array
        for(int i=0; i<aF1star.length; i++){
            aF1star[i].E(((IntegratorVelocityVerlet.MyAgent)atomAgentMin.getAgent((IAtomLeaf)listMin.getAtom(i))).force());
            aF2star[i].Ea1Tv1(2.0, aF[i]);
            aF2star[i].ME(aF1star[i]);
        }
    }
	
	/**
	 * Finds the component of the dimer's total force that is perpendicular to the normal vector.
	 * 
	 * @param aN Array of vectors that create an n-dimensional vector parallel to the orientation of the dimer.
	 * @param aF1 Array of vectors holding the forces in box min
	 * @param aFperp Array of vectors that create an n-dimensional vector perpendicular to the orientation of the dimer.
	 */
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

	public Object makeAgent(IAtomLeaf a) {
		return new IntegratorVelocityVerlet.MyAgent(space);
	}

	public void releaseAgent(Object agent, IAtomLeaf atom) {
		// TODO Auto-generated method stub	
	}
	
}
