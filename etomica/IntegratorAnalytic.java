package etomica;

/**
 * Integrator that generates atom trajectories from an analytic formula.
 * Takes an IntegratorAnalytic.AtomAction instance and performs action on all atoms in each
 * time step; intention is for the action to set the position of the atom for the
 * current time (but it could do anything).
 *
 * @author David Kofke
 * @author Nancy Cribbin
 */
 
 /* History of changes
  * 7/31/02 (DAK) new
  */
 
 public class IntegratorAnalytic extends IntegratorMD {
    
    private final AtomIteratorList atomIterator = new AtomIteratorList();
    private AtomAction action;
    
    public IntegratorAnalytic() {
        this(Simulation.instance);
    }
    public IntegratorAnalytic(Simulation sim) {
        super(sim);
    }
    
    public void doStep() {
        if(action == null) return;
        action.setTime(elapsedTime());
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom atom = atomIterator.next();
            action.actionPerformed(atom);
        }
    }
    
    public void setAction(AtomAction action) {this.action = action;}
    
    public AtomAction getAction() {return action;}
    
	/**
	 * Overrides superclass method to instantiate iterators when iteratorFactory in phase is changed.
	 * Called by Integrator.addPhase and Integrator.iteratorFactorObserver.
	 */
	public boolean addPhase(Phase p) {
	    if(!super.addPhase(p)) return false;
        atomIterator.setBasis(p.speciesMaster.atomList);
        return true;
    }
    
    public Integrator.Agent makeAgent(Atom a) {return null;}
    
    
    /**
     * Extends AtomAction class to add a method to set the time.
     */
    public static abstract class AtomAction extends etomica.AtomAction {
        protected double time;
        public void setTime(double t) {time = t;}
    }
    
 }//end of IntegratorAnalytic
 