package etomica;
import etomica.units.Dimension;
import etomica.units.Unit;
import etomica.units.Count;

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
                                                Integrator.IntervalListener.BeforePbc, 
                                                Integrator.IntervalListener.AfterPbc, 
                                                EtomicaElement {

    public static final String getVersion() {return "MeterMeanSquareDisplacementFixed:02.09.08/"+MeterScalar.VERSION;}
 
    private int nAtoms = 0;
    AtomIterator iterator;
    private Space.Vector[][] rDelta;
    private Space.Vector[] rLast;
    private Space.Vector rAccum;
    private int iPtr;
    private double timeStep;
    
    public MeterMeanSquareDisplacementFixed() {
        this(Simulation.instance);
    }
    public MeterMeanSquareDisplacementFixed(AtomIterator iter) {
        this(Simulation.instance);
        setAtoms(iter);
    }
    public MeterMeanSquareDisplacementFixed(Simulation sim) {
        super(sim);
        setLabel("Mean square displacement");
        setUpdateInterval(1);
        rAccum = sim.space.makeVector();
        setActive(true);
        setXLabel("time");
        setXMin(0.0);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Mean squared displacement, suitable for diffusion calculations");
        return info;
    }
    
    public etomica.units.Dimension getXDimension() {return etomica.units.Dimension.TIME;}

    public Unit defaultIOUnit() {return new Unit(Count.UNIT);}
    
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
        rDelta = new Space.Vector[nAtoms][nPoints];
        rLast = new Space.Vector[nAtoms];
        iter.reset();
        int i=0;
        while(iter.hasNext()) {
            for(int j=0; j<nPoints; j++) {
                rDelta[i][j] = simulation().space.makeVector();
            }
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
    
    public double[] currentValue() {
            
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
    /**
     *  Override superclass method to perform update of atom displacements.
     */
    public void intervalAction(Integrator.IntervalEvent evt) {
        if(!active) return;
        if(evt.type() != Integrator.IntervalEvent.INTERVAL) return; //don't act on start, done, initialize events
	    //ignore iieCount; want to update every interval
	    
        iterator.reset();
        int i = 0;
        //evaluate difference from last coordinate before pbc applied
        if(evt.isBeforePbc()) { 
            while(iterator.hasNext()) {
                Space.Vector r = iterator.next().coord.position();
                rDelta[i][iPtr].Ev1Mv2(r,rLast[i]);
                i++;
            }
        }
        //store last coordinate after pbc applied
        else { 
           while(iterator.hasNext()) {rLast[i++].E(iterator.next().coord.position());}
           super.intervalAction(evt);
        }   
    }//end of intervalAction	    
}//end of class