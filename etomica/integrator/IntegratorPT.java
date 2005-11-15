package etomica.integrator;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataSource;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveListener;
import etomica.integrator.mcmove.MCMoveSwapConfiguration;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
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

public class IntegratorPT extends IntegratorManagerMC implements EtomicaElement {
    
    public IntegratorPT(PotentialMaster potentialMaster) {
        this(potentialMaster, MCMoveSwapConfiguration.FACTORY);
    }
    
    public IntegratorPT(PotentialMaster potentialMaster, MCMoveSwapFactory swapFactory) {
        super(potentialMaster);
        setGlobalMoveInterval(100);
        mcMoveSwapFactory = swapFactory;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Parallel-tempering Monte Carlo simulation");
        return info;
    }
    
    /**
     * Adds the given integrator to those managed by this integrator, and
     * includes integrator's phase to the set among which configurations are
     * swapped.  Phase of new integrator will be subject to swapping with
     * the most recently added integrator/phase (and the next one, if another
     * is added subsequently).
     */
	public void addIntegrator(IntegratorPhase integrator){
	    super.addIntegrator(integrator);
        
		if (nIntegrators > 1) {
            MCMove newMCMove = mcMoveSwapFactory.makeMCMoveSwap(potential, (IntegratorPhase)integrators[nIntegrators-1], 
                    (IntegratorPhase)integrators[nIntegrators]);
            moveManager.addMCMove(newMCMove);
		}
	}
    
	private final MCMoveSwapFactory mcMoveSwapFactory;

	
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
	    public MCMove makeMCMoveSwap(PotentialMaster potentialMaster, IntegratorPhase integrator1, IntegratorPhase integrator2);
	}

    /**
     * Interface for a move that swaps two phases.  Enables access to
     * the swapped phases.
     */
    public interface MCMoveSwap {
        public Phase[] swappedPhases();
    }


	/**
     * Meter that tracks the swapping of the phases in parallel-tempering
     * simulation.  Designed for input to a DisplayPlot to provide a graphical
     * record of how the phases swap configurations.
     */
    public static class PhaseTracker implements DataSource, MCMoveListener, java.io.Serializable {
        
        private int[] track;
        private double[] dtrack;
        private DataDoubleArray data;
        
        public PhaseTracker() {
            data = new DataDoubleArray("Phase Tracker", Dimension.NULL,0);
        }
        
        public DataInfo getDataInfo() {
            return data.getDataInfo();
        }
        
        /**
         * Method called when two phases are successfully exchanged.
         */
        public void actionPerformed(MCMoveEvent evt) {
            if(evt.isTrialNotify || !evt.wasAccepted) return;
            if(!(evt.mcMove instanceof MCMoveSwap)) return;
            Phase[] phases = ((MCMoveSwap)evt.mcMove).swappedPhases();
            int i0 = phases[0].getIndex()-1;
            int i1 = phases[1].getIndex()-1;
            int temp = track[i0];
            track[i0] = track[i1];
            track[i1] = temp;
        }
        
        /**
         * Specifies the number of phases that are tracked.
         */
        public void setNumPhases(int numPhases) {
            track = new int[numPhases];
            data = new DataDoubleArray(data.getDataInfo(), new int[] {numPhases});
            dtrack = data.getData();
            for(int i=0; i<numPhases; i++) {
                track[i] = i;
            }
        }
        
        /**
         * Returns array y such that y[i] is the current
         * phase of the configuration that began in phase i.
         */
        public Data getData() {
            for (int i=0; i<track.length; i++) {
                dtrack[i] = track[i];
            }
            return data;
        }
        
        public int getDataLength() {
            return track.length;
        }
        
    }
    
}

