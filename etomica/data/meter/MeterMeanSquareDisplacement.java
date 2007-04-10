package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorPhase;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.units.Undefined;

/**
 *  Computes the mean square displacement for a set of atoms.
 *  Use of this meter usually involves another meter must calling this meter's
 *  currentValue method, which then reports the current mean square displacement.
 *  This meter does not perform any averaging of its own.
 *
 *  @author David Kofke
 */

//seriously consider using MeterMeanSquareDisplacementFixed instead (in development project)

//doesn't implement Meter because phase information comes from integrator
public class MeterMeanSquareDisplacement extends DataSourceScalar {

    
    public MeterMeanSquareDisplacement(Space space, IntegratorPhase integrator) {
        this(space, integrator, new AtomIteratorLeafAtoms());
    }

    public MeterMeanSquareDisplacement(Space space, IntegratorPhase integrator, AtomIteratorPhaseDependent iter) {
        super("Mean square displacement", Undefined.DIMENSION);
        this.space = space;
        this.integrator = integrator;
        setIterator(iter);
        integrator.addListener(new BeforePbc(this));
        integrator.addListener(new AfterPbc(this));
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Mean squared displacement, suitable for diffusion calculations");
        return info;
    }
    
    public void setIterator(AtomIteratorPhaseDependent iterator) {
        if(iterator == null) {
            throw new NullPointerException("Cannot give a null iterator");
        }
        this.iterator = iterator;
        if(integrator.getPhase() != null) {
            iterator.setPhase(integrator.getPhase());
            reset();
        } else {
            //throw an exception, because meter won't be informed when integrator has phase set
            throw new IllegalStateException("Must first define Phase for Integrator before constructing MeterMeanSquareDisplacement");
        }
    }

    /**
     * Specifies the set of atoms that will be tracked to compute their MSD.
     * The given iterator loops through the atoms.
     */
    public void reset() {
        nAtoms = iterator.size();
        rAccum = new IVector[nAtoms];
        rLast = new IVector[nAtoms];
        iterator.reset();
        int i=0;
        while(iterator.hasNext()) {
            rAccum[i] = space.makeVector();
            rLast[i] = space.makeVector();
            rLast[i].E(((AtomLeaf)iterator.nextAtom()).getCoord().getPosition());
            i++;
        }
    }
    
    public double getDataAsScalar() {
        double sum = 0.0;
        for(int i=0; i<nAtoms; i++) {
            sum += rAccum[i].squared();
        }
        return sum/nAtoms;
    }
    
    private static class BeforePbc implements IntegratorIntervalListener, java.io.Serializable {
        BeforePbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public int getPriority() {return 50;}//PBC is 100-199
        public void intervalAction() {
//            meter.iterator.setPhase(meter.integrator.getPhase()[0]);
            meter.iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            while(meter.iterator.hasNext()) {
                IVector r = ((AtomLeaf)meter.iterator.nextAtom()).getCoord().getPosition();
                meter.rAccum[i].PE(r);
                meter.rAccum[i].ME(meter.rLast[i]);
                meter.rLast[i].E(r);
                i++;
            }
        }
        private static final long serialVersionUID = 1L;
        final MeterMeanSquareDisplacement meter;
    }
    
    private static class AfterPbc implements IntegratorIntervalListener, java.io.Serializable {
        AfterPbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public int getPriority() {return 200;}//PBC is 100-199
        public void intervalAction() {
            meter.iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            //store last coordinate after pbc applied
           while(meter.iterator.hasNext()) {meter.rLast[i++].E(((AtomLeaf)meter.iterator.nextAtom()).getCoord().getPosition());}
        }
        private static final long serialVersionUID = 1L;
        final MeterMeanSquareDisplacement meter;
    }
    
    private static final long serialVersionUID = 1L;
    private int nAtoms = 0;
    AtomIteratorPhaseDependent iterator;
    IntegratorPhase integrator;
    protected IVector[] rAccum, rLast;
    private final Space space;

}//end of class
