package etomica.data.meter;
import etomica.AtomIterator;
import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.Space;
import etomica.integrator.IntegratorMD;
import etomica.space.Vector;
import etomica.units.Count;
import etomica.units.Dimension;
import etomica.units.Unit;

/**
 *  Computes the mean square displacement for a set of atoms for a fixed
 *  time interval.  Keeps track of positions of targeted atoms for a period 
 *  of time, and evaluates MSD for all time differences within this period.
 *  Differs from MeterMeanSquareDisplacement, which keeps MSD for a (typically
 *  larger) set of atoms from a starting point to the current time; then 
 *  MSD vs time is given by History.  In the present class, MSD vs time
 *  is the function kept by this MeterFunction.
 *
 *  @author David Kofke
 */
 
 /* History
  * 09/08/02 (DAK) new from MeterMeanSquareDisplacement.
  */

public class MeterMeanSquareDisplacementFixed extends MeterFunction implements 
                                                EtomicaElement {

    public static final String getVersion() {return "MeterMeanSquareDisplacementFixed:02.09.08/"+etomica.data.meter.VERSION;}
 
    private int nAtoms = 0;
    AtomIterator iterator;
    private Vector[][] rDelta;
    private Vector[] rLast;
    private Vector rAccum;
    private int iPtr;
    private double timeStep;
    
    public MeterMeanSquareDisplacementFixed(Space space, AtomIterator iter) {
        this(space);
        setAtoms(iter);
    }
    public MeterMeanSquareDisplacementFixed(Space space) {
        super();
        setLabel("Mean square displacement");
        setUpdateInterval(1);
        rAccum = space.makeVector();
        setActive(true);
        setXLabel("time");
        setXMin(0.0);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Mean squared displacement, suitable for diffusion calculations");
        return info;
    }
    
    public etomica.units.Dimension getXDimension() {return etomica.units.Dimension.TIME;}

    public Unit defaultIOUnit() {return Count.UNIT;}
    
    public void setPhaseIntegrator(Integrator integrator) {
        super.setPhaseIntegrator(integrator);
        timeStep = ((IntegratorMD)integrator).getTimeStep();
        setXMax(nPoints*timeStep);
    }
    
    /**
     * Override to check that integrator's time step has not changed.
     * If so, x values must be updated first.
     */
    public double[] xValues() {
        if(integrator != null && ((IntegratorMD)integrator).getTimeStep() != timeStep) {
            timeStep = ((IntegratorMD)integrator).getTimeStep();
            setXMax(nPoints*timeStep);
        }
        return super.xValues();
    }
    
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
        iPtr = 0;
        nAtoms = iter.size();
        rDelta = new Vector[nAtoms][nPoints];
        rLast = new Vector[nAtoms];
        iter.reset();
        int i=0;
        while(iter.hasNext()) {
            for(int j=0; j<nPoints; j++) {
                rDelta[i][j] = simulation().space.makeVector();
            }
            rLast[i] = simulation().getSpace().makeVector();
            rLast[i].E(iter.next().coord.position());
            i++;
        }
    }
    public AtomIterator getAtoms() {return iterator;}
    
    public void reset() {
        super.reset();
        if(iterator != null) setAtoms(iterator);
    }
    
    public double[] getData() {
            
        for(int j=0; j<nPoints; j++) {
            y[j] = 0.0;
        }
        for(int i=0; i<nAtoms; i++) {
            rAccum.E(0.0);
            int idx = iPtr;
            for(int j=0; j<nPoints; j++) {
                rAccum.PE(rDelta[i][idx--]);
                y[j] += rAccum.squared();
                if(idx < 0) idx = nPoints-1;
            }
        }
        for(int j=0; j<nPoints; j++) {
            y[j] /= (double)nAtoms;
        }
        
        //advance position where next deltaR is stored
        iPtr++;
        if(iPtr == nPoints) iPtr = 0;
        
        return y;
    }
    
    private class BeforePbc implements Integrator.IntervalListener {
        public int getPriority() {return 50;}//PBC is 100-199
        public void intervalAction(Integrator.IntervalEvent evt) {
            if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
            iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            while(iterator.hasNext()) {
                Vector r = iterator.nextAtom().coord.position();
                rDelta[i][iPtr].Ev1Mv2(r,rLast[i]);
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