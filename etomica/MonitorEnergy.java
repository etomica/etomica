package etomica;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveListener;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.space3d.Space3D;
import etomica.units.Dimension;

/**
 * Tracks the value of the energy of the configuration.  Updates
 * value with each Monte Carlo trial, so that the total energy
 * of the system can be known without requiring it be recomputed
 * from scratch.
 *
 * @author David Kofke
 */
 
 /* History of changes
  *  7/10/02 new
  */
public class MonitorEnergy implements MCMoveListener, PhaseListener, java.util.Observer, java.io.Serializable {
    
    public MonitorEnergy() {
        this(Simulation.instance);
    }
    
    public MonitorEnergy(Simulation sim) {
        setRefreshInterval(10000);
        setWarningTolerance(1e-2);
        setIncludeLrc(true);
        energy = sim.energySum(this);
        potential = sim.hamiltonian.potential;
    }
    
    /**
     * Sets the phase for which the energy is monitored.  Causes this
     * to be a listener for reset() calls in the phase, registers it with the
     * integrator monitor in the phase (so that this is informed if the phase's
     * integrator is changed) and registers this as a MCMove listener of 
     * the phase's integrator, if it has one.
     */
    public void setPhase(Phase phase) {
        if(this.phase != null) this.phase.integratorMonitor.deleteObserver(this);
        this.phase = phase;
        phase.integratorMonitor.addObserver(this);
        phase.addListener(this);
        if(phase.integrator() != null) setIntegrator((IntegratorMC)phase.integrator());
    }
    
    /**
     * Called when integrator is set in phase.  Causes update of integrator.
     * Implementation of Observer interface.
     */
    public void update(java.util.Observable integratorMonitor, Object integrator) {
        setIntegrator((IntegratorMC)integrator);
    }
    
    private void setIntegrator(IntegratorMC integrator) {
        if(!(integrator instanceof IntegratorMC)) 
            throw new IllegalArgumentException("Attempt to apply MonitorEnergy to phase with integrator that is not IntegratorMC");
        if(integrator != null) integrator.removeMCMoveListener(this);
        this.integrator = (IntegratorMC)integrator;
        integrator.addMCMoveListener(this);
    }
    
    /**
     * Checks if MCMove has generated a new configuration, and if so adds
     * incremental change in energy to current energy.
     */
    public void actionPerformed(MCMoveEvent evt) {
        if(evt.isTrialNotify) return;
        if(evt.wasAccepted) currentValue += evt.mcMove.energyChange(phase);
        if(--counter == 0) refresh();
    }
    
    /**
     * Recomputes the current value of the energy from scratch in the event of
     * a call to reset in the monitored phase.  Implementation of PhaseListener 
     * interface.
     */
    public void actionPerformed(PhaseEvent evt) {
        if(evt.type == PhaseEvent.RESET) {
            currentValue = potential.calculate(phase, iteratorDirective, energy.reset()).sum();
        }
    }
    
    /**
     * Recalculates energy, and compares with current value; issues warning 
     * (writes to console) if fractional difference is greater than warningTolerance.
     */
    public void refresh() {
        counter = refreshInterval;
        double oldValue = currentValue;
        currentValue = potential.calculate(phase, iteratorDirective, energy.reset()).sum();
        if(Math.abs((oldValue-currentValue)/oldValue) > warningTolerance) {
            System.out.println("Warning in MonitorEnergy.refresh");
            System.out.println("Old, new energy: "+oldValue+"  "+currentValue);
        }
        System.out.println("Old, new energy, difference: "+oldValue+"  "+currentValue+"  "+(currentValue-oldValue));
    }
    
    /**
     * Sets the number of energy-update increments that are performed before
     * a call to refresh() to recalculate the energy from scratch.  Default
     * value is 10000.
     */
    public void setRefreshInterval(int i) {
        if(i <= 0) return;
        refreshInterval = i;
        counter = refreshInterval;
    }
    /**
     * Accessor method for the refresh interval.
     */
    public int getRefreshInterval() {return refreshInterval;}
    /**
     * Indicates that refresh interval has dimensions of QUANTITY.
     */
    public Dimension getRefreshIntervalDimension() {return Dimension.QUANTITY;}
    
    /**
     * Sets the threshold value of the error in the current and refreshed
     * energy; if the fractional error is above this threshold, a warning
     * is issued in refresh.
     */
    public void setWarningTolerance(double eps) {warningTolerance = eps;}
    /**
     * Accessor method for warning tolerance.
     */
    public double getWarningTolerance() {return warningTolerance;}
    /**
     * Indicates that the warning tolerance is dimensionless.
     */
    public Dimension getWarningToleranceDimension() {return Dimension.NULL;}
    
    /**
     * Mutator method for flag that determines if long-range correction to
     * potential should be included in energy.  Default is <code>true</code>
     */
    public void setIncludeLrc(boolean b) {iteratorDirective.includeLrc = b;}
    /**
     * Accessor method for flag indicating whether long-range correction is
     * included in energy.
     */
    public boolean isIncludeLrc() {return iteratorDirective.includeLrc;}
    
    /**
     * Returns the (stored) value of the energy of the current configuration
     * of the phase being monitored by this instance.  Returns NaN until a
     * phase is specified (via setPhase) and a configuration is set up in it.
     */
    public double currentValue() {return currentValue;}
    
    private double currentValue = Double.NaN;
    private Phase phase;
    private IntegratorMC integrator;
    private int counter, refreshInterval;
    private double warningTolerance;
    private final IteratorDirective iteratorDirective = new IteratorDirective();
    private final PotentialCalculationEnergySum energy;
    private final PotentialMaster potential;
    
    /**
     * Method to demonstrate and test the use of this class.
     */
    public static void main(String[] args) {
        Simulation sim = new Simulation(new Space3D());
      //  setIteratorFactory(new IteratorFactoryCell(this));
        Simulation.instance = sim;
	    Phase phase = new Phase(sim);
	    Potential2 potential = new P2LennardJones();
	    Species species = new SpeciesSpheresMono(sim);
	    potential.setSpecies(species, species);
	    Controller controller = new Controller(sim);
	    IntegratorMC integrator = new IntegratorMC(sim);
	    new MCMoveAtom(integrator);
	    new MCMoveVolume(integrator);
	    MonitorEnergy  monitor = new MonitorEnergy(sim);
	    monitor.setPhase(phase);
	    monitor.setRefreshInterval(1000);
//	    sim.elementCoordinator.go();
	    controller.start();
	    
    }    
}//end of MonitorEnergy