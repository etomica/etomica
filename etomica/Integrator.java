package etomica;

import java.util.EventObject;
import java.util.Vector;

import etomica.units.Dimension;
import etomica.utility.NameMaker;

/**
 * Integrator is used to define the algorithm used to move the atoms around and
 * generate new configurations in one or more phases. All integrator methods,
 * such as molecular dynamics or Monte Carlo are implemented via subclasses of
 * this Integrator class. The Integrator's activities are managed via the
 * actions of the governing Controller.
 * 
 * @author David Kofke
 */

/*
 * History
 * 
 * 07/10/03 (DAK) made Agent interface public 08/25/03 (DAK) changed default for
 * doSleep to <false> 01/27/04 (DAK) initialized iieCount to inverval (instead
 * of interval+1) in run method; changed setInterval do disallow non-positive
 * interval 04/13/04 (DAK) modified reset such that doReset is called if running
 * is false
 */
public abstract class Integrator implements java.io.Serializable {

	protected final PotentialMaster potential;

	protected Phase firstPhase;

	protected Phase[] phase;

	protected boolean equilibrating = false;

	int phaseCount = 0;

	int phaseCountMax = 1;

	protected int sleepPeriod = 10;

	//should use a different collection structure
	private Vector intervalListenersBeforePbc = new Vector();

	private Vector intervalListenersImposePbc = new Vector();

	private Vector intervalListenersAfterPbc = new Vector();

	int integrationCount = 0;

	protected double temperature = Default.TEMPERATURE;

	protected boolean isothermal = false;
    
    private String name;

	public Integrator(PotentialMaster potentialMaster) {
        setName(NameMaker.makeName(this.getClass()));
		phase = new Phase[phaseCountMax];
		this.potential = potentialMaster;
        if (Default.AUTO_REGISTER) {
            Simulation.getDefault().register(this);
        }
	}


    /**
     * Accessor method of the name of this phase
     * 
     * @return The given name of this phase
     */
    public final String getName() {return name;}
    /**
     * Method to set the name of this simulation element. The element's name
     * provides a convenient way to label output data that is associated with
     * it.  This method might be used, for example, to place a heading on a
     * column of data. Default name is the base class followed by the integer
     * index of this element.
     * 
     * @param name The name string to be associated with this element
     */
    public void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the phase
     */
    public String toString() {return getName();}

    /**
	 * Performs the elementary integration step, such as a molecular dynamics
	 * time step, or a Monte Carlo trial.
	 */
	public abstract void doStep();

	/**
	 * Defines the actions taken by the integrator to reset itself, such as
	 * required if a perturbation is applied to the simulated phase (e.g.,
	 * addition or deletion of a molecule). Also invoked when the
	 * integrator is started or initialized.
	 */
	public abstract void reset(); 
	
	/**
	 * Returns a new instance of an agent of this integrator for placement in
	 * the given atom in the ia (IntegratorAgent) field.
	 */
	public abstract Object makeAgent(Atom a);

	/**
	 * Initializes the integrator, performing the following steps: (1) deploys
	 * agents in all atoms; (2) call doReset method; (3) fires an event
	 * indicating to registered listeners indicating that initialization has
	 * been performed (i.e. fires IntervalEvent of type field set to
	 * INITIALIZE).
	 */
	public void initialize() {
		deployAgents();
		reset();
	}

	//how do agents get placed in atoms made during the simulation?
	protected void deployAgents() { //puts an Agent of this integrator in each
									// atom of all phases
		AtomIteratorListSimple iterator = new AtomIteratorListSimple();
		for (int i = 0; i < phaseCount; i++) {
			Phase p = phase[i];
			iterator.setList(p.speciesMaster.atomList);
			iterator.reset();
			while (iterator.hasNext()) {//does only leaf atoms; do atom groups
										// need agents?
				Atom a = iterator.nextAtom();
				a.setIntegratorAgent(makeAgent(a));
			}
		}
	}

	public void setTemperature(double t) {
		temperature = t;
	}

	public final double getTemperature() {
		return temperature;
	}

	public final double temperature() {
		return temperature;
	}

	public Dimension getTemperatureDimension() {
		return Dimension.TEMPERATURE;
	}

	//Other introspected properties
	public void setIsothermal(boolean b) {
		isothermal = b;
	}

	public boolean isIsothermal() {
		return isothermal;
	}

	/**
	 * @return Returns flag indicating whether integrator is in equilibration mode.
	 */
	public boolean isEquilibrating() {
		return equilibrating;
	}

	/**
	 * @param equilibrating
	 *            Sets equilibration mode of integrator.
	 */
	public void setEquilibrating(boolean equilibrating) {
		this.equilibrating = equilibrating;
	}

	/**
	 * @return true if integrator can perform integration of another phase,
	 *         false if the integrator has all the phases it was built to handle
	 */
	public boolean wantsPhase() {
		return phaseCount < phaseCountMax;
	}

	public Phase getPhase(int i) {
		return phase[i];
	}

	/**
	 * Performs activities needed to set up integrator to work on given phase.
	 * This method should not be called directly; instead it is invoked by the
	 * phase in its setIntegrator method.
	 * 
	 * @return true if the phase was successfully added to the integrator; false
	 *         otherwise
	 */
	//perhaps should throw an exception rather than returning a boolean "false"
	public boolean addPhase(Phase p) {
		for (int i = 0; i < phaseCount; i++) {
			if (phase[i] == p)
				return false;
		} //check that phase is not already registered
		if (!this.wantsPhase()) {
			return false;
		} //if another phase not wanted, return false
		phase[phaseCount] = p;
		phaseCount++;
		firstPhase = phase[0];
		return true;
	}

	/**
	 * Performs activities needed to disconnect integrator from given phase.
	 * This method should not be called directly; instead it is invoked by the
	 * phase in its setIntegrator method
	 */
	public void removePhase(Phase p) {
		for (int i = 0; i < phaseCount; i++) {
			if (phase[i] == p) {//phase found; remove it
				phase[i] = null;
				phaseCount--;
				if (phaseCount > 0)
					phase[i] = phase[phaseCount];
				firstPhase = phase[0];
				break;
			}
		}
	}

	/**
	 * Attempts to set the given integrator as the integrator for all phases
	 * that were added to this phase.
	 */
	public void transferPhasesTo(Integrator anotherIntegrator) {
		for (int i = 0; i < phaseCount; i++) {
			phase[i].setIntegrator(anotherIntegrator);
		}
	}

	public synchronized void addIntervalListener(IntervalListener iil) {
		boolean added = false;
		//must check all possibilities because listener may implement multiple
		// pbc interfaces
		if (iil instanceof IntervalListener.BeforePbc) {
			intervalListenersBeforePbc.addElement(iil);
			added = true;
		}
		if (iil instanceof IntervalListener.ImposePbc) {
			intervalListenersImposePbc.addElement(iil);
			added = true;
		}
		if (iil instanceof IntervalListener.AfterPbc) {
			intervalListenersAfterPbc.addElement(iil);
			added = true;
		}
		//if not implementing any of the pbc interfaces, default is afterPbc
		if (!added)
			intervalListenersAfterPbc.addElement(iil);
	}

	public synchronized void removeIntervalListener(IntervalListener iil) {
		boolean removed = false;
		if (iil instanceof IntervalListener.BeforePbc) {
			intervalListenersBeforePbc.removeElement(iil);
			removed = true;
		}
		if (iil instanceof IntervalListener.ImposePbc) {
			intervalListenersImposePbc.removeElement(iil);
			removed = true;
		}
		if (iil instanceof IntervalListener.AfterPbc) {
			intervalListenersAfterPbc.removeElement(iil);
			removed = true;
		}
		if (!removed)
			intervalListenersAfterPbc.removeElement(iil);
	}

	/**
	 * Notifies registered listeners that an interval has passed. Not
	 * synchronized, so unpredictable behavior if listeners are added while
	 * notification is in process (this should be rare).
	 */
	public void fireIntervalEvent(IntervalEvent iie) {
		iie.setBeforePbc(true);
		int n = intervalListenersBeforePbc.size();
		for (int i = 0; i < n; i++) {
			IntervalListener listener = (IntervalListener) intervalListenersBeforePbc
					.elementAt(i);
			listener.intervalAction(iie);
		}
		n = intervalListenersImposePbc.size();
		for (int i = 0; i < n; i++) {
			IntervalListener listener = (IntervalListener) intervalListenersImposePbc
					.elementAt(i);
			listener.intervalAction(iie);
		}
		iie.setBeforePbc(false);
		n = intervalListenersAfterPbc.size();
		for (int i = 0; i < n; i++) {
			IntervalListener listener = (IntervalListener) intervalListenersAfterPbc
					.elementAt(i);
			listener.intervalAction(iie);
		}
	}

	/**
	 * Registers with the given integrator all listeners currently registered
	 * with this integrator. Removes all listeners from this integrator.
	 */
	public synchronized void transferListenersTo(Integrator anotherIntegrator) {
		if (anotherIntegrator == this)
			return;
		int n = intervalListenersBeforePbc.size();
		for (int i = 0; i < n; i++) {
			IntervalListener listener = (IntervalListener) intervalListenersBeforePbc
					.elementAt(i);
			anotherIntegrator.addIntervalListener(listener);
		}
		n = intervalListenersImposePbc.size();
		for (int i = 0; i < n; i++) {
			IntervalListener listener = (IntervalListener) intervalListenersImposePbc
					.elementAt(i);
			anotherIntegrator.addIntervalListener(listener);
		}
		n = intervalListenersAfterPbc.size();
		for (int i = 0; i < n; i++) {
			IntervalListener listener = (IntervalListener) intervalListenersAfterPbc
					.elementAt(i);
			anotherIntegrator.addIntervalListener(listener);
		}
		intervalListenersBeforePbc.removeAllElements();
		intervalListenersImposePbc.removeAllElements();
		intervalListenersAfterPbc.removeAllElements();
	}

	public Phase[] getPhases() {
		return phase;
	}

	/**
	 * Integrator agent that holds a force vector. Used to indicate that an atom
	 * could be under the influence of a force.
	 */
	public interface Forcible {
		public Space.Vector force();
	}

	public static class IntervalEvent extends EventObject {

		//Typed constants used to indicate the type of event integrator is
		// announcing
		public static final Type START = new Type("Start"); //simulation is
															// starting

		public static final Type INTERVAL = new Type("Interval"); //routine
																  // interval
																  // event

		public static final Type DONE = new Type("Done"); //simulation is
														  // finished

		public static final Type INITIALIZE = new Type("Initialize"); //integrator
																	  // is
																	  // initializing

		private final Type type;

		private boolean beforePbc;

		private int interval;

		public IntervalEvent(Integrator source, Type t) {
			super(source);
			type = t;
		}

		public IntervalEvent(Integrator source, int interval) {
			this(source, INTERVAL);
			this.interval = interval;
		}

		public int getInterval() {
			return interval;
		}

		/**
		 * Indicates if notification is before or after the PhaseImposePbc
		 * listeners have been notified. Returns true if notifying BeforePbc or
		 * PhaseImposePbc listeners; return false if notifying AfterPbc
		 * listeners.
		 */
		public final boolean isBeforePbc() {
			return beforePbc;
		}

		/**
		 * Sets the before/after status relative to imposePbc listeners. Should
		 * be used only by Integrator that is firing event.
		 */
		final void setBeforePbc(boolean b) {
			beforePbc = b;
		}

		public Type type() {
			return type;
		}

		//class used to mark the different types of interval events
		private final static class Type extends Constants.TypedConstant {
			private Type(String label) {
				super(label);
			}

			public static final Constants.TypedConstant[] choices = new Constants.TypedConstant[] {
					START, INTERVAL, DONE, INITIALIZE };

			public final Constants.TypedConstant[] choices() {
				return choices;
			}
		}
	}

	public interface IntervalListener extends java.util.EventListener {
		public void intervalAction(IntervalEvent evt);

		/**
		 * Marker interface that indicates an IntervalListener that should be
		 * notified before any application of periodic boundary conditions.
		 */
		public interface BeforePbc extends IntervalListener {
		}

		/**
		 * Marker interface that indicates an IntervalListener that will invoke
		 * periodic boundary conditions when notified.
		 */
		public interface ImposePbc extends IntervalListener {
		}

		/**
		 * Marker interface that indicates an IntervalListener that should be
		 * notified after any application of periodic boundary conditions. This
		 * is the default for any IntervalListener that doesn't have a
		 * IntervalListenerPBC marker interface.
		 */
		public interface AfterPbc extends IntervalListener {
		}
	}

}

