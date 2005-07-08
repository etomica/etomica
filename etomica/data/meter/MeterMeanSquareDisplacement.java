package etomica.data.meter;
import etomica.AtomIterator;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;
import etomica.Meter;
import etomica.Space;
import etomica.space.Vector;
import etomica.units.Dimension;

/**
 *  Computes the mean square displacement for a set of atoms.
 *  Use of this meter usually involves another meter must calling this meter's
 *  currentValue method, which then reports the current mean square displacement.
 *  This meter does not perform any averaging of its own.
 *
 *  @author David Kofke
 */

public class MeterMeanSquareDisplacement extends Meter implements 
                                                EtomicaElement {

    private int nAtoms = 0;
    AtomIterator iterator;
    private Vector[] rAccum, rLast;
    private Vector deltaR;
    private final Space space;
    
    public MeterMeanSquareDisplacement(Space space, AtomIterator iter) {
        this(space);
        setAtoms(iter);
    }
    public MeterMeanSquareDisplacement(Space space) {
        super(1);
        this.space = space;
        setLabel("Mean square displacement");
        deltaR = space.makeVector();
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Mean squared displacement, suitable for diffusion calculations");
        return info;
    }

//    public Unit defaultIOUnit() {return Count.UNIT;}
    
    /**
     * Returns dimensions of this meter's output (current UNDEFINED).
     */
    public Dimension getDimension() {return Dimension.UNDEFINED;}
    
    /**
     * Specifies the set of atoms that will be tracked to compute their MSD.
     * The given iterator loops through the atoms.
     */
    public void setAtoms(AtomIterator iter) {
        iterator = iter;
        nAtoms = iter.size();
        rAccum = new Vector[nAtoms];
        rLast = new Vector[nAtoms];
        iter.reset();
        int i=0;
        while(iter.hasNext()) {
            rAccum[i] = space.makeVector();
            rLast[i] = space.makeVector();
            rLast[i].E(iter.nextAtom().coord.position());
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
    
    private class BeforePbc implements IntegratorIntervalListener, java.io.Serializable {
        public int getPriority() {return 50;}//PBC is 100-199
        public void intervalAction(IntegratorIntervalEvent evt) {
             iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            while(iterator.hasNext()) {
                Vector r = iterator.nextAtom().coord.position();
                rAccum[i].PE(r);
                rAccum[i].ME(rLast[i]);
                rLast[i].E(r);
                i++;
            }
        }//end of intervalAction    
    }//end of BeforePbc
    
    private class AfterPbc implements IntegratorIntervalListener, java.io.Serializable {
        public int getPriority() {return 200;}//PBC is 100-199
        public void intervalAction(IntegratorIntervalEvent evt) {
            iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            //store last coordinate after pbc applied
           while(iterator.hasNext()) {rLast[i++].E(iterator.nextAtom().coord.position());}
        }//end of intervalAction    
    }//end of AfterPbc
    
}//end of class
