package etomica;

import java.awt.Graphics;
import java.util.Observable;
import java.util.Observer;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.utility.Histogram;

/**
 * Superclass of all meter types.
 * A meter is responsible for measuring and averaging properties during a simulation.
 * A meter is placed in a phase with the setPhase method, and it does is measurements exclusively on that phase,
 * although there are exceptions.<br>
 * Meters do not carry any ability to display their measurements.  Display of meter data is accomplished
 * by connecting a meter to an appropriate Display object.<br>
 * All meter measurements are performed in simulation units.  Conversion to other units is completed
 * by the Display object that is presenting the meter data.<br>
 * The abstract meter is subclassed into a Meter and a MeterFunction.  An extension of the Meter subclass 
 * measures a scalar property (such as the energy), while an extension of the MeterFunction subclass
 * measures a vector or function property (such as the radial distribution function).<br>
 * A meter updates its averages after receiving some number of interval events from an integrator.
 * All of the manipulations involving the calculation of averages and error bars are handled by
 * the Accumulator inner class of MeterAbstract.  Each meter contains one (or more, for a MeterFunction)
 * instance of an accumulator.<br>
 * A meter may be set to keep a histogram of its measured values, in addition to the average and variance.
 * 
 *
 * @author David Kofke
 * @see Meter
 * @see MeterFunction
 * @see MeterMultiFunction
 * @see MeterAbstract.Accumulator
 */
public abstract class MeterAbstract implements Integrator.IntervalListener, Simulation.Element,  java.io.Serializable
{
    /**
     * Number of integration interval events received before another call to updateSums
     */
    protected int updateInterval;
    /**
     * Counter that keeps track of the number of interval events received since last call to updateSums
     */
    protected int iieCount;
    /**
     * The phase in which this meter is performing its measurements
     */
    protected Phase phase;
    /**
     * The integrator firing the IntervalEvents telling the meter to take its measurements
     */
     protected Integrator integrator;
    /**
     * A string describing the property measured by the meter
     */
    protected String label;
    /**
     * Flag specifying whether the meter responds to integrator events
     * If false, the meter does not perform regular measurements or keep averages
     * In this situation the meter is probably measuring a property for use by some other object
     * Default is <code>true</code> for a Meter, but <code>false</code> for a MeterFunction.
     */
    protected boolean active;
    
    private String name;
    private final Simulation parentSimulation;
    private boolean added = false;
    
    private transient Observer integratorObserver;
    private transient Observer boundaryObserver;
    private transient Observer iteratorFactoryObserver;
    
    private boolean histogramming = false;

	public MeterAbstract(Simulation sim) {
	    parentSimulation = sim;
	    setUpdateInterval(1);
	    label = "Property";
	    parentSimulation.register(this);
	}
	
    public final Simulation parentSimulation() {return parentSimulation;}
    public final Class baseClass() {return MeterAbstract.class;}
    public final boolean wasAdded() {return added;}
    public final void setAdded(boolean b) {added = b;}
    	
	/**
	 * Returns the physical dimensions (e.g., mass, length, pressure, etc.) of the quantity being measured
	 */
	public abstract Dimension getDimension();
	
	/**
	 * Returns the default unit of measurement for the metered quantity.
	 * Obtains this from the default unit prescribed by the dimension for this quantity,
	 * which in turn is given by the Simulation.unitSystem field.
	 */
	public Unit defaultIOUnit() {return getDimension().defaultIOUnit();}
	
	/**
	 * Method to set reset (discard accumulated values for) all averages and other sums
	 */
	public abstract void reset();
	/**
	 * Method called upon receiving updateInterval interval-events from the integrator
	 */
	public abstract void updateSums();

    /**
     * Method to render an image of the meter to the screen.  
     * Not commonly used, and default is to perform no action
     */
    public void draw(Graphics g, int[] origin, double scale) {} 
    /**
     * Sets the phase in which the meter is performing its measurements.
     * If given phase is null, meter is removed from current phase, if already in one.
     * Meter is reset if new phase is not null.
     */
	public void setPhase(Phase p) {
	    //remove observers from current phase, if not null
	    if(phase != null) {
//	        phase.removeMeter(this);
    	    if(integratorObserver != null) phase.integratorMonitor.deleteObserver(integratorObserver);
    	    if(boundaryObserver != null) phase.boundaryMonitor.deleteObserver(boundaryObserver);
	        if(iteratorFactoryObserver != null) phase.iteratorFactoryMonitor.deleteObserver(iteratorFactoryObserver);
	    }
	    if(p == null) { //setting meter to have no phase
	        setPhaseIntegrator(null);
	        phase = null; //this must be done after setPhaseIntegrator call
	        return;
	    }
	    else { //given phase is not null
	        //update the phase handle
	        phase = p;
//	        p.addMeter(this);
	        //list this meter as an interval listener of phase integrator
	        setPhaseIntegrator(p.integrator()); 
	        //add an observer to the phase to make sure this meter is notified if the integrator changes
	        integratorObserver = integratorObserver();
	        p.integratorMonitor.addObserver(integratorObserver);
	        //perform similar steps to register and observe phase boundary and iteratorFactory if so indicated
	        if(usesPhaseBoundary()) {
    	        boundaryObserver = boundaryObserver();
    	        setPhaseBoundary(p.boundary());
	            p.boundaryMonitor.addObserver(boundaryObserver);
	        }
	        if(usesPhaseIteratorFactory()) {
	            iteratorFactoryObserver = iteratorFactoryObserver();
    	        setPhaseIteratorFactory(p.iteratorFactory());
	            p.iteratorFactoryMonitor.addObserver(iteratorFactoryObserver);
	        }
	        //zero all averages for this meter
	        reset();
	    }
	}
	
	/**
	 * Sets the integrator for this meter.
	 * Registers the meter as an interval listener to the integrator, and if
	 * appropriate, also as a collision listener.
	 * This method is invoked by the setPhase method and, if the phase changes its
	 * integrator, by an observer that is registered within the phase as part of the setPhase process.
	 */
	protected void setPhaseIntegrator(Integrator newIntegrator) {
	    //note that phase does not update its integrator until after notifying observers, 
	    //so this should work to remove this meter from interval-listeners list of old integrator
	    Integrator oldIntegrator = null;
	    if(phase != null) oldIntegrator = phase.integrator();
	    if(oldIntegrator != null) {
	        oldIntegrator.removeIntervalListener(this);
	        if(oldIntegrator instanceof IntegratorHard && this instanceof IntegratorHard.CollisionListener) 
	            ((IntegratorHard)oldIntegrator).removeCollisionListener((IntegratorHard.CollisionListener)this);
	    }
	    if(newIntegrator != null) {
	        newIntegrator.addIntervalListener(this);
	        if(newIntegrator instanceof IntegratorHard && this instanceof IntegratorHard.CollisionListener) 
	            ((IntegratorHard)newIntegrator).addCollisionListener((IntegratorHard.CollisionListener)this);
	    }
	    this.integrator = newIntegrator;
	}
	    
    /**
     * The integrator firing the interval events that tell this meter to do its update.
     * Integrator is set by the setPhaseIntegrator method, which is called by setPhase.
     */
     public Integrator integrator() {return integrator;}

	/**
	 * Returns an observer that can be registered with the phase's integratorMonitor.
	 * In this way this meter is informed if the integrator object changes to another instance in the phase.
	 * Registration of this observer is done automatically by the setPhase method.
	 */
	protected Observer integratorObserver() {
	    return new Observer() {
	        //This is the action that is to happen if phase takes a new integrator
	        public void update(Observable o, Object arg) {
	            setPhaseIntegrator((Integrator)arg);  //if looking for a bug, check that this calls overridden method in subclasses
	        }
	    };
	}
	/**
	 * Returns an observer that can be registered with the phase's iteratorFactoryMonitor.
	 * In this way this meter is informed if the iteratorFactory object changes to another instance in the phase.
	 * Registration of this observer is done automatically by the setPhase method, but 
	 * only if the usesIteratorFactory method of this meter returns <code>true</code>.
	 */
	protected Observer iteratorFactoryObserver() {
	    return new Observer() {
	        //This is the action that is to happen if phase takes a new iteratorFactory
	        public void update(Observable o, Object arg) {
	            setPhaseIteratorFactory((IteratorFactory)arg);
	        }
	    };
	}
	/**
	 * Returns an observer that can be registered with the phase's boundaryMonitor
	 * In this way this meter is informed if the boundary object changes to another instance in the phase
	 * Registration of this observer is done automatically by the setPhase method, but 
	 * only if the usesBoundary method of this meter returns <code>true</code>.
	 */
	protected Observer boundaryObserver() {
	    return new Observer() {
	        //This is the action that is to happen if phase takes a new boundary
	        public void update(Observable o, Object arg) {
	            setPhaseBoundary((Space.Boundary)arg);
	        }
	    };
	}
	/**
	 * Flags whether the meter uses the phase boundary when making its measurements.
	 * If so, then the meter will be registered (within the setPhase method) as an observer of the phase's boundaryMonitor,
	 * which notifies if the boundary object in the phase is changed to another boundary.
	 * This method is declared abstract, rather than given a default value, to minimize 
	 * typographical errors in subclasses that try to declare it to return <code>true</code>.
	 */
	protected abstract boolean usesPhaseBoundary();
	/**
	 * Flags whether the meter uses the phase iteratorFactory when making its measurements.
	 * If so, then the meter will be registered (within the setPhae method) as an observer of the phase's iteratorFactoryMonitor,
	 * which notifies if the iteratorFactory object in the phase is changed to another iteratorFactory.
	 * This method is declared abstract, rather than given a default value, to minimize 
	 * typographical errors in subclasses that try to declare it to return <code>true</code>.
	 */
	protected abstract boolean usesPhaseIteratorFactory();
	
	/**
	 * Method performs no action, but can be overridden in subclasses to handle setting or change of boundary in phase
	 */
	protected void setPhaseBoundary(Space.Boundary b) {}
	/**
	 * Method performs no action, but can be overridden in subclasses to handle setting or change of iteratorFactory in phase
	 */
	protected void setPhaseIteratorFactory(IteratorFactory i) {}
	
	/**
	 * Accessor method for the phase in which this meter resides
	 */
	 public Phase getPhase() {return phase;}
	
	/**
	 * Accessor method for the meter's label
	 */
	public String getLabel() {return label;}
	/**
	 * Accessor method for the meter's label
	 */
	public void setLabel(String s) {label = s;}

    /**
     * Accessor method for the updateInterval
     * @see #updateInterval
     */
    public final int getUpdateInterval() {return updateInterval;}
    /**
     * Accessor method for the updateInterval
     * Sets to given value and resets count of interval events
     * @see #updateInterval
     */
    public final void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
    
    /**
     * Sets whether the meter responds to integrator events
     * If false, the meter does not perform regular measurements or keep averages
     */
    public void setActive(boolean b) {active = b;}
    /**
     * Returns status of active flag
     * @see setActive
     */
    public final boolean isActive() {return active;}
    
    /**
     * Integrator.IntervalListenter interface method.
     * Counts number of events received and calls updateSums after receiving event updateInterval times
     */
    public void intervalAction(Integrator.IntervalEvent evt) {
        if(!active) return;
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        updateSums();
	    }
    }
    
	/**
	* Accessor method to indicate if the meter should keep a histogram of all measured values.
	* Default is false (do not keep histogram).
	*/
	public abstract boolean isHistogramming();
    	 
	/**
	* Accessor method to indicate if the meter should keep a histogram of all measured values.
	* Default is false (do not keep histogram).
	*/
	public abstract void setHistogramming(boolean b);

    /**
     * Accessor method of the name of this object
     * 
     * @return The given name
     */
    public final String getName() {return name;}

    /**
     * Method to set the name of this object
     * 
     * @param name The name string to be associated with this object
     */
    public final void setName(String name) {this.name = name;}

    /**
     * Overrides the Object class toString method to have it return the output of getName
     * 
     * @return The name given to the object
     */
    public String toString() {return getName();}  //override Object method
          
    /**
     * Class for accumulation of averages and error statistics
     * One or more instances of this class is present in every meter.
     * Presently error value is not computed well (assumes independent measurements)
     * Future development project: include block averaging, add a plot display to show error as a function of block size
     */
    public final static class Accumulator implements java.io.Serializable {
        private double sum, sumSquare;
        private double mostRecent = Double.NaN;
        private int count;
        private Histogram histogram;
        private boolean histogramming = false;
        
        public Accumulator() {
            reset();
        }
        
        /**
         * Add the given value to the sum and sum of squares
         */
        public void add(double value) {
            mostRecent = value;
            if(Double.isNaN(value)) return;
	        sum += value;
	        sumSquare += value*value;
	        if(histogramming) histogram.addValue(value);
	        count++;
        }
        
        /**
         * Compute and return the average from the present value of the accumulation sum
         */
	    public double average() {
	        return (count>0) ? sum/(double)count : Double.NaN;
	    }
    	
    	/**
    	 * Computer and return the 67% confidence limits of the average
    	 *   Needs more development
    	 */
	    public double error() {    //temporary---needs to be rewritten to do block averaging
	        double avg = average();
	        return (count>1) ? Math.sqrt((sumSquare/(double)count - avg*avg)/(double)(count-1)) : Double.NaN;
	    }
	    
	    /**
	     * Compute and return the variance of the recorded data.
	     */
	     public double variance() {
	        double avg = sum/(double)count;
	        return (count>0) ? sumSquare/(double)count - avg*avg : Double.NaN;
	     }
	    
	    /**
	     * Returns the value last passed to the add method
	     */
	    public double mostRecent() {return mostRecent;}
    	
    	/**
    	 * Resets all sums to zero
    	 */
	    public void reset() {
	        count = 0;
	        sum = 0;
	        sumSquare = 0;
	        if(histogram != null) histogram.reset();
	    }
	    /**
	     * Returns the histogram in its current state.
	     */
	    public Histogram histogram() {return histogram;}
	    
	    /**
	    * Accessor method to indicate if the meter should keep a histogram of all measured values.
	    * Default is false (do not keep histogram).
	    */
	    public boolean isHistogramming() {return histogramming;}
    	 
	    /**
	    * Accessor method to indicate if the meter should keep a histogram of all measured values.
	    * Default is false (do not keep histogram).
	    */
	    public void setHistogramming(boolean b) {
	        histogramming = b;
	        if(histogramming && histogram == null) histogram = new Histogram();
	    }
    	 
	}
}	 
