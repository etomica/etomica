package etomica;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.units.Count;

/**
 *  Computes the mean square displacement for a set of atoms.
 *  Use of this meter usually involves another meter must calling this meter's
 *  currentValue method, which then reports the current mean square displacement.
 *  This meter does not perform any averaging of its own.
 *
 *  @author David Kofke
 */

public class MeterMeanSquareDisplacement extends MeterAbstract implements 
                                                Integrator.IntervalListener.BeforePbc, 
                                                Integrator.IntervalListener.AfterPbc, 
                                                EtomicaElement {

    public static final String getVersion() {return "MeterMeanSquareDisplacement:01.03.24/"+MeterScalar.VERSION;}
 
    private int nAtoms = 0;
    AtomIterator iterator;
    private Space.Vector[] rAccum, rLast;
    private Space.Vector deltaR;
    
    public MeterMeanSquareDisplacement(Space space, AtomIterator iter) {
        this(space);
        setAtoms(iter);
    }
    public MeterMeanSquareDisplacement(Space space) {
        super(1);
        setLabel("Mean square displacement");
        setUpdateInterval(1);
        deltaR = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Mean squared displacement, suitable for diffusion calculations");
        return info;
    }

    public Unit defaultIOUnit() {return Count.UNIT;}
    
    /**
     * Overrides superclass method so that updateInterval is always 1.
     * Meter will not work properly if it skips updating at each interval event,
     * because PBC could be invoked by another IntervalListener, and this would
     * invalidate mean square displacement calculation.
     */
    public void setUpdateInterval(int i) {super.setUpdateInterval(1);}
        
    /**
     * Returns dimensions of this meter's output, which in this case is QUANTITY.
     */
    public Dimension getDimension() {return Dimension.NULL;}
    
    /**
     * Specifies the set of atoms that will be tracked to compute their MSD.
     * The given iterator loops through the atoms.
     */
    public void setAtoms(AtomIterator iter) {
        iterator = iter;
        nAtoms = iter.size();
        rAccum = new Space.Vector[nAtoms];
        rLast = new Space.Vector[nAtoms];
        iter.reset();
        int i=0;
        while(iter.hasNext()) {
            rAccum[i] = simulation().space().makeVector();
            rLast[i] = simulation().space().makeVector();
            rLast[i].E(iter.next().coord.position());
            i++;
        }
    }
    public AtomIterator getAtoms() {return iterator;}
    
    public void reset() {
        super.reset();
        if(iterator != null) setAtoms(iterator);
    }
    
    public void doMeasurement() {
        double sum = 0.0;
        for(int i=0; i<nAtoms; i++) {
            sum += rAccum[i].squared();
        }
        data[0] = sum/(double)nAtoms;
    }
    /**
     *  Override superclass method to perform update of atom displacements.
     *  Update sums is not called, since simple averaging of mean square displacement
     *  is not of interest.  
     */
    public void intervalAction(Integrator.IntervalEvent evt) {
        if(!active) return;
        if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
	    //ignore iieCount; want to update every interval
	    
        iterator.reset();
        int i = 0;
        //accumulate difference from last coordinate before pbc applied
        if(evt.isBeforePbc()) { 
            while(iterator.hasNext()) {
                Space.Vector r = iterator.next().coord.position();
                rAccum[i].PE(r);
                rAccum[i].ME(rLast[i]);
                rLast[i].E(r);
                i++;
            }
        }
        //store last coordinate after pbc applied
        else { 
           while(iterator.hasNext()) {rLast[i++].E(iterator.next().coord.position());}
        }   
	    if(--iieCount == 0) {
	        iieCount = updateInterval;
	        accumulator.history().addValue(getData());
	    }
    }//end of intervalAction	    
}//end of class