package etomica;

import etomica.units.Dimension;

/**
 * Parallel-tempering integrator.  Oversees other integrators that are defined to perform
 * MC trials (or perhaps molecular dynamics) in different phases.  These integrators
 * are identified (added) to this integrator, and the simulation runs on this
 * integrator's thread.  When this integrator does a step, it passes on the
 * instruction to each added integrator, causing all to do a step on each's phase.
 * Occasionally, this integrator will instead attempt to swap the configurations
 * of two of the phases.  Acceptance of this move depends on parameters (usually 
 * temperature) of the integrators for the phases, as well as the configurations
 * in each phase.  The swap trial is performed by a MCMove class designed for
 * this purpose.  Such a class is made by a MCMoveSwapConfigurationFactory class
 * that is identified to this integrator (a default is selected if not specified).
 * Every time an integrator is added to this one, a MCMoveSwap class is made (by this
 * integrator using the factory) to manage swap trials between the new integrator's
 * phase and that of the one most recently added.
 * 
 * @author David Kofke
 */
 
 /* History of changes
  * 7/16/02 (DAK) new
  */

public class IntegratorPT extends IntegratorMC implements EtomicaElement {
    
    public IntegratorPT() {
        this(Simulation.instance);
    }
    public IntegratorPT(SimulationElement parent) {
        this(parent, new MCMoveSwapFactoryDefault());
    }
    public IntegratorPT(SimulationElement parent, MCMoveSwapFactory swapFactory) {
        super(parent);
        setSwapInterval(100);
        mcMoveSwapFactory = swapFactory;
    }
    
    /**
     * Adds the given integrator to those managed by this integrator, and
     * includes integrator's phase to the set among which configurations are
     * swapped.  Phase of new integrator will be subject to swapping with
     * the most recently added integrator/phase (and the next one, if another
     * is added subsequently).
     */
	public void addIntegrator(Integrator integrator){
		int n = integrators.length;
		Integrator[] newintegrators = new Integrator[n+1];
		for(int i=0; i<n; i++) {newintegrators[i] = integrators[i];}
		newintegrators[n] = integrator;
		integrators = newintegrators;
		nIntegrators++;
		
		if(nIntegrators > 1) {
		    MCMove[] newMCMoveSwap = new MCMove[n];
		    for(int i=0; i<n-1; i++) {newMCMoveSwap[i] = mcMoveSwap[i];}
		    newMCMoveSwap[n-1] = mcMoveSwapFactory.makeMCMoveSwap(this, integrators[n-1], integrators[n]);
            mcMoveSwap = newMCMoveSwap;
		}
	}
	
	/**
	 * Fires interval event for this integrator, then instructs
	 * each sub-integrator to fire event.
	 */
	public void fireIntervalEvent(Integrator.IntervalEvent ie) {
	    super.fireIntervalEvent(ie);
	    for(int i=0; i<nIntegrators; i++) {
	        integrators[i].fireIntervalEvent(ie);
	    }
	}
	
	
	/**
     * Performs a Monte Carlo trial that attempts to swap the configurations
     * between two "adjacent" phases, or instructs all integrators to perform
     * a single doStep.
     */
    public void doStep() {
        if(Simulation.random.nextDouble() < swapProbability) {
            super.doStep();
        } else {
            for(int i=0; i<nIntegrators; i++) integrators[i].doStep();
        }
    }
    
    /**
     * Resets this integrator and passes on the reset to all managed integrators.
     */
    public void reset() {
        super.reset();
	    for(int i=0; i<nIntegrators; i++) {
	        integrators[i].reset();
	    }
	}
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Parallel-tempering Monte Carlo simulation");
        return info;
    }
    
    /**
     * Returns an array containing all the MCMove classes that perform
     * swaps between the phases.
     */
    public MCMove[] swapMoves() {return mcMoveSwap;}
    
    /**
     * Sets the average interval between phase-swap trials.  With each 
     * call to doStep of this integrator, there will be a probability of
     * 1/nSwap that a swap trial will be attempted.  A swap is attempted
     * for only one pair of phases.  Default is 100.
     */
    public void setSwapInterval(int nSwap) {
        if(nSwap < 1) nSwap = 1;
        this.nSwap = nSwap;
        swapProbability = 1.0/(double)nSwap;
    }
    
    /**
     * Accessor method for the average interval between phase-swap trials.
     */
    public int getSwapInterval() {return nSwap;}
    
    private int nSwap;
    private double swapProbability;
	private Integrator[] integrators = new Integrator[0];
	private MCMove[] mcMoveSwap = new MCMove[0];
	private int nIntegrators = 0;
	private final MCMoveSwapFactory mcMoveSwapFactory;

	
	//----------- inner interface -----------
	
	/**
	 * Interface for a class that can make a MCMove that will swap
	 * the configurations of two phases.  Different MCMove classes
	 * would do this differently, depending on ensemble of simulation
	 * and other factors.
	 */
public interface MCMoveSwapFactory {
    
    /**
     * @param integratorMC the parent integrator using this move
     * @param integrator1 integrator for one of the phases being swapped
     * @param integrator2 integrator for the other phase
     */
    public MCMove makeMCMoveSwap(IntegratorMC integratorMC, Integrator integrator1, Integrator integrator2);
}//end of MCMoveSwapFactory

    // --------------- inner interface ------------
    
    /**
     * Interface for a move that swaps two phases.  Enables access to
     * the swapped phases.
     */
public interface MCMoveSwap {
    
    public Phase[] swappedPhases();
    
}



    // -----------inner class----------------
    
public static class MCMoveSwapFactoryDefault implements MCMoveSwapFactory {
    public MCMove makeMCMoveSwap(IntegratorMC integratorMC, 
                                    Integrator integrator1, Integrator integrator2) {
        return new MCMoveSwapConfiguration(integratorMC, integrator1, integrator2);
    }
}//end of MCMoveSwapFactoryDefault 


	
	// -----------inner class----------------
	
	/**
	 * Basic MCMove for swapping coordinates of atoms in two phases.
	 * Requires same number of atoms in each phase.
	 */
public static class MCMoveSwapConfiguration extends MCMove implements MCMoveSwap {

	private Integrator integrator1, integrator2;	
    private final IteratorDirective iteratorDirective = new IteratorDirective();
	private AtomIteratorAllMolecules iterator1 = new AtomIteratorAllMolecules();
	private AtomIteratorAllMolecules iterator2 = new AtomIteratorAllMolecules();
	private AtomIteratorAllMolecules affectedAtomIterator = new AtomIteratorAllMolecules();
	private Space.Vector r;
	private double u1, u2, temp1, temp2, deltaU1;
	private Phase phase1, phase2;
	private final Phase[] swappedPhases = new Phase[2];

	public MCMoveSwapConfiguration(IntegratorMC integratorMC, 
	                                Integrator integrator1, Integrator integrator2) {
  		super(integratorMC);
		r = integratorMC.space.makeVector();
		setTunable(false);
		this.integrator1 = integrator1;
		this.integrator2 = integrator2;
	}

	public boolean doTrial() {
 		phase1 = integrator1.getPhase(0);
		phase2 = integrator2.getPhase(0);

		temp1 = integrator1.getTemperature();
		temp2 = integrator2.getTemperature();

        u1 = potential.calculate(phase1, iteratorDirective, energy.reset()).sum();
        u2 = potential.calculate(phase2, iteratorDirective, energy.reset()).sum();
        deltaU1 = Double.NaN;
        return true;
    }
    
    public double lnTrialRatio() {return 0.0;}
    
    public double lnProbabilityRatio() {
        deltaU1 = u2 - u1;  //if accepted, energy of phase1 will be u2, and its old energy is u1
		return  -deltaU1*((1/temp1) - (1/temp2));
	}
	
	/**
	 * Swaps positions of molecules in two phases.
	 */
	public void acceptNotify() {
		iterator1.setBasis(phase1);
		iterator2.setBasis(phase2);
			
		iterator1.reset();
		iterator2.reset();

		while(iterator1.hasNext()) {
			Atom a1 = iterator1.next();
			Atom a2 = iterator2.next();

			r.E(a1.coord.position());
				
			a1.coord.translateTo(a2.coord.position());
			a2.coord.translateTo(r);
		}
	}
	
	/**
     * Performs no action; nothing required when move rejected.
     */
	public void rejectNotify() {}
	
	/**
	 * Implementation of MCMoveSwap interface
	 */
	public Phase[] swappedPhases() {
	    swappedPhases[0] = phase1;
	    swappedPhases[1] = phase2;
	    return swappedPhases;
	}

	public double energyChange(Phase phase) {
	    if(phase == phase1)      return +deltaU1;
	    else if(phase == phase2) return -deltaU1;
	    else                     return 0.0;
	}
	
	public AtomIterator affectedAtoms(Phase p) {
	    if(p == phase1 || p == phase2) {
	        affectedAtomIterator.setBasis(p);
	        affectedAtomIterator.reset();
	        return affectedAtomIterator;
	    }
	    else return AtomIterator.NULL;
	}
}//end of MCMoveSwapConfiguration

//----------------- inner class ------------------------------

/**
 * Meter that tracks the swapping of the phases in parallel-tempering
 * simulation.  Designed for input to a DisplayPlot to provide a graphical
 * record of how the phases swap configurations.
 */
//XXX this does not seem to fit the criteria for a Meter.  It is not associated with a phase
    public class MeterPhaseTracker extends MeterFunction implements MCMoveListener {
        
        private DataSource[] histories;
        private int[] track;
        
        public MeterPhaseTracker() {
            super(IntegratorPT.this.simulation());
            setActive(true);
            IntegratorPT.this.addMCMoveListener(this);
        }
        
        /**
         * Method called when two phases are successfully exchanged.
         */
        public void actionPerformed(MCMoveEvent evt) {
            if(evt.isTrialNotify || !evt.wasAccepted) return;
            if(!(evt.mcMove instanceof MCMoveSwap)) return;
            Phase[] phases = ((MCMoveSwap)evt.mcMove).swappedPhases();
            int i0 = phases[0].index;
            int i1 = phases[1].index;
            int temp = track[i0];
            track[i0] = track[i1];
            track[i1] = temp;
        }
        
        public Dimension getXDimension() {return Dimension.NULL;}
        public Dimension getDimension() {return Dimension.NULL;}
        
        /**
         * Specifies phases that are tracked.
         */
         //don't need to pass phases; just number of phases
        public void setPhases(Phase[] trackedPhases) {
            setNPoints(trackedPhases.length);
            track = new int[trackedPhases.length];
            setHistorying(true);
            histories = new DataSource[trackedPhases.length];
            for(int i=0; i<nPoints; i++) {
                track[i] = i;
                histories[i] = accumulator[i].history();
            }
        }
        
        /**
         * Returns data source array that is suitable to input to
         * setDataSource method of DisplayPlot.
         */
        public DataSource[] dataSource() {return histories;}
        
        /**
         * Returns array y such that y[i] is the current
         * phase of the configuration that began in phase i.
         */
        public double[] getData() {
            for(int i=0; i<track.length; i++) {
                y[track[i]] = i;
            }
            return y;
        }
        
    }//end of MeterPhaseTracker
    
}//end of IntegratorPT

