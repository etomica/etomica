package etomica;

/**
 * Meter for the root-mean-square velocity of a set of atoms.  
 * Useful to obtain histograms of atom speeds.
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 08/30/02 (DAK) new
  */
  
public class MeterVelocityRms extends MeterAbstract implements MeterAtomic {
    
    public MeterVelocityRms() {
        this(Simulation.instance);
    }
    public MeterVelocityRms(Simulation sim) {
        super(sim, 1);
    }
    
    /**
     * Returns the average of the current value over all atoms.
     */
    public void doMeasurement() {
        iterator.reset();
        double sum = 0.0;
        while(iterator.hasNext()) {
            sum += currentValue(iterator.next());
        }
        data[0] = sum/iterator.size();
    }
    
    /**
     * Not implemented.
     */
    public etomica.units.Dimension getDimension() {
        return etomica.units.Dimension.UNDEFINED;
    }
    
    /**
     * Returns the speed of the given atom.
     */
    public double currentValue(Atom a) {
        return a.coord.rm()*Math.sqrt(a.coord.momentum().squared());
    }
    
    /**
     * Overrides superclass method to select all leaf atoms in phase for iteration,
     * but only if iterator was not previously set (or was set to null).
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
        if(iterator == AtomIterator.NULL) iterator = p.makeAtomIterator();
    }
    
    /** 
     * Overrides superclass method to add data looping over all atoms.
     */
	public void updateSums() {
	    iterator.reset();
	    while(iterator.hasNext()) {
	        accumulator.add(currentValue(iterator.next()));
	    }
	}
	
	public void setAtoms(AtomIterator iterator) {
	    if(iterator == null) this.iterator = AtomIterator.NULL;
	    else this.iterator = iterator;
	}
	public AtomIterator getAtoms() {return iterator;}

    private AtomIterator iterator = AtomIterator.NULL;
    
}//end of MeterVelocityRms