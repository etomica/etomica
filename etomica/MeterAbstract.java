package etomica;

import java.awt.Graphics;
import java.util.Observable;
import java.util.Observer;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.utility.Histogram;
import etomica.utility.History;

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
public abstract class MeterAbstract implements Integrator.IntervalListener, Simulation.Element,  java.io.Serializable {

    public static final String VERSION = "MeterAbstract:01.05.14";
    
    /**
     * Number of integration interval events received before another call to updateSums
     */
    protected int updateInterval;
    /**
     * Number of interval events received between additions to history.
     */
   // protected int historyInterval;
    /**
     * Counter that keeps track of the number of interval events received since last call to updateSums
     */
    protected int iieCount;
    /**
     * Counter that keeps track of the interval events since last history accumulation.
     */
   // protected int historyCount;
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
    
    protected int historyWindow = 100;
    
    /**
     * Size of subaveraging block used to evaluate confidence limits.
     * Default is 1000 updateIntervals.
     */
    int blockSize = Default.BLOCK_SIZE;
    
    private String name;
    private final Simulation parentSimulation;
    private boolean added = false;
    
    private transient Observer integratorObserver;
    private transient Observer boundaryObserver;
    private transient Observer iteratorFactoryObserver;
    
    boolean histogramming = false;
    boolean historying = false;

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
	 * Method called upon receiving updateInterval interval-events from the integrator
	 */
	public abstract void updateSums();

	protected abstract Accumulator[] allAccumulators();
 
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
    public void setUpdateInterval(int i) {
        if(i > 0) {
            updateInterval = i;
            iieCount = updateInterval;
        }
    }
    
    /**
     * BlockSize describes the size of the subaverages used to estimate confidence limits 
     * on simulation averages.
     * Default is 1000 updateIntervals.
     */
    public void setBlockSize(int b) {blockSize = Math.max(1,b);}
    /**
     * BlockSize describes the size of the subaverages used to estimate confidence limits 
     * on simulation averages.
     */
    public int getBlockSize() {return blockSize;}
    
    
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
        if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        updateSums();
	    }
    }
    
	/**
	 * Method to set reset (discard accumulated values for) all averages and other sums
	 */
	public void reset() {
	    Accumulator[] accumulators = allAccumulators();
	    for(int i=0; i<accumulators.length; i++) {
	        accumulators[i].reset();
	    }
	}

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
	    Accumulator[] accumulators = allAccumulators();
	    for(int i=0; i<accumulators.length; i++) {
	        accumulators[i].setHistogramming();
	    }
	}

	/**
	* Accessor method to indicate if the meter should keep a history of all measured values.
	* Default is false (do not keep history).
	*/
	public boolean isHistorying() {return historying;}
    	 
	/**
	* Accessor method to indicate if the meter should keep a history of all measured values.
	* Default is false (do not keep history).
	*/
	public void setHistorying(boolean b) {
	    historying = b;
	    Accumulator[] accumulators = allAccumulators();
	    for(int i=0; i<accumulators.length; i++) {
	        accumulators[i].setHistorying();
	    }
	}
	
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
    public final void setName(String name) {
        this.name = name;
	    Accumulator[] accumulators = allAccumulators();
	    if(accumulators == null) return;
	    for(int i=0; i<accumulators.length; i++) {
	        if(accumulators[i].histogram() != null) 
	            accumulators[i].histogram().setName(this.toString() + ":Histogram");
	        if(accumulators[i].history() != null) 
	            accumulators[i].history().setName(this.toString() + ":History");
	    }
    }

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
    public final class Accumulator implements java.io.Serializable {
        private double sum, sumSquare, blockSum, error;
        private double mostRecent = Double.NaN;
        private double mostRecentBlock = Double.NaN;
        private int count, blockCountDown;
        private Histogram histogram;
        private History history;
        
        public Accumulator() {
            reset();
        }
        
        /**
         * Add the given value to the sum and block sum
         */
        public void add(double value) {
            mostRecent = value;
            if(Double.isNaN(value)) return;
            
	        blockSum += value;
	        if(--blockCountDown == 0) {//count down to zero to determine completion of block
	            blockSum /= blockSize;//compute block average
	            sum += blockSum;
	            sumSquare += blockSum*blockSum;
    	        count++;
    	        if(count > 1) {
    	            double avg = sum/(double)count;
    	            error = Math.sqrt((sumSquare/(double)count - avg*avg)/(double)(count-1));
    	        }
	            //reset blocks
	            mostRecentBlock = blockSum;
	            blockCountDown = blockSize;
	            blockSum = 0.0;
	        }
	        if(histogramming) histogram.addValue(value);
	        if(historying) history.addValue(value);
        }
        
        /**
         * Compute and return the average from the present value of the accumulation sum
         */
	    public double average() {
	        int blockCount = blockSize - blockCountDown;
	        return (count+blockCount > 0) ? 
	            (sum + blockSum/blockSize)/(double)(count + (double)blockCount/(double)blockSize) 
	            : Double.NaN;
	    }
    	
    	/**
    	 * Return the 67% confidence limits of the average based on variance in block averages.
    	 */
	    public double error() {
	        return error;
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
    	 * Returns the value of the most recent block average.
    	 */
    	public double mostRecentBlock() {
    	    return (mostRecentBlock!=Double.NaN) ? mostRecentBlock : blockSum/(blockSize - blockCountDown);
        }
        
    	/**
    	 * Resets all sums to zero
    	 */
	    public void reset() {
	        count = 0;
	        sum = 0.0;
	        sumSquare = 0.0;
	        error = Double.NaN;
	        blockCountDown = blockSize;
	        blockSum = 0.0;
	        if(histogram != null) histogram.reset();
	        if(history != null) history.reset();
	    }
	    
	    /**
	     * Returns the histogram in its current state.
	     */
	    public Histogram histogram() {return histogram;}
	    
	    /**
	     * Returns the history of values recorded by the accumulator.
	     */
	    public History history() {return history;}
    	 
	    /**
	    * Accessor method to indicate if the meter should keep a histogram of all measured values.
	    * Default is false (do not keep histogram).
	    * Does not take an argument; true/falue value is obtained from outer-class meter
	    */
	    //need way to update name if meter name changes
	    public void setHistogramming() {
	        if(histogramming && histogram == null) {
	            histogram = new Histogram();
	            histogram.setName(MeterAbstract.this.toString() + ":Histogram");
	            histogram.setLabel(MeterAbstract.this.getLabel() + " histogram");
	            histogram.setXDimension(MeterAbstract.this.getDimension());
	        }
	    }
	    
	    /**
	    * Accessor method to indicate if the meter should keep a history of all measured values.
	    * Default is false (do not keep history).
	    * Does not take an argument; true/falue value is obtained from outer-class meter
	    */
	    public void setHistorying() {
	        if(historying && history == null) {
	            history = new History();
	            history.setName(MeterAbstract.this.toString() + ":History");
	            history.setLabel(MeterAbstract.this.getLabel() + " history");
	        }
	    }
	    
	}//end of Accumulator
	
	/**
	 * Typed constant that can be used to indicated the quantity
	 * to be taken from a meter (e.g., average, error, current value, etc.).
	 * Used primarily by Display objects.
	 */
	public static class ValueType extends DataSource.ValueType {
        public ValueType(String label) {super(label);}
        
        public Constants.TypedConstant[] choices() {return CHOICES;}
        public static final ValueType AVERAGE = new ValueType("Average");
        public static final ValueType ERROR = new ValueType("67% Confidence Limits");
        public static final ValueType CURRENT = new ValueType("Current value");
        public static final ValueType MOST_RECENT = new ValueType("Latest value");
        public static final ValueType MOST_RECENT_BLOCK = new ValueType("Latest block average");
        public static final ValueType VARIANCE = new ValueType("Variance");
        public static final Constants.TypedConstant[] CHOICES = 
            new ValueType[] {
                AVERAGE, ERROR, MOST_RECENT, MOST_RECENT_BLOCK, VARIANCE};
	}
	
}//end of MeterAbstract	 
