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
                                                EtomicaElement {

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
    
    private class BeforePbc implements Integrator.IntervalListener {
        public int getPriority() {return 50;}//PBC is 100-199
        public void intervalAction(Integrator.IntervalEvent evt) {
            if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
            iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            while(iterator.hasNext()) {
                Space.Vector r = iterator.nextAtom().coord.position();
                rAccum[i].PE(r);
                rAccum[i].ME(rLast[i]);
                rLast[i].E(r);
                i++;
            }
        }//end of intervalAction    
    }//end of BeforePbc
    
    private class AfterPbc implements Integrator.IntervalListener {
        public int getPriority() {return 200;}//PBC is 100-199
        public void intervalAction(Integrator.IntervalEvent evt) {
            if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
            iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            //store last coordinate after pbc applied
           while(iterator.hasNext()) {rLast[i++].E(iterator.nextAtom().coord.position());}
        }//end of intervalAction    
    }//end of AfterPbc
    
}//end of class