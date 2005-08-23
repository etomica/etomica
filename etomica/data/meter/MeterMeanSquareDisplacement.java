package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.Integrator;
import etomica.IntegratorIntervalEvent;
import etomica.IntegratorIntervalListener;
import etomica.Space;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceScalar;
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

//seriously consider using MeterMeanSquareDisplacementFixed instead (in development project)

//doesn't implement Meter because phase information comes from integrator
public class MeterMeanSquareDisplacement extends DataSourceScalar {

    
    public MeterMeanSquareDisplacement(Space space, Integrator integrator) {
        this(space, integrator, new AtomIteratorLeafAtoms());
    }

    public MeterMeanSquareDisplacement(Space space, Integrator integrator, AtomIteratorPhaseDependent iter) {
        super("Mean square displacement", Dimension.UNDEFINED);
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
        if(integrator.getPhase().length > 0 && integrator.getPhase()[0] != null) {
            iterator.setPhase(integrator.getPhase()[0]);
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
        rAccum = new Vector[nAtoms];
        rLast = new Vector[nAtoms];
        iterator.reset();
        int i=0;
        while(iterator.hasNext()) {
            rAccum[i] = space.makeVector();
            rLast[i] = space.makeVector();
            rLast[i].E(iterator.nextAtom().coord.position());
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
        public void intervalAction(IntegratorIntervalEvent evt) {
//            meter.iterator.setPhase(meter.integrator.getPhase()[0]);
            meter.iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            while(meter.iterator.hasNext()) {
                Vector r = meter.iterator.nextAtom().coord.position();
                meter.rAccum[i].PE(r);
                meter.rAccum[i].ME(meter.rLast[i]);
                meter.rLast[i].E(r);
                i++;
            }
        }//end of intervalAction    
        final MeterMeanSquareDisplacement meter;
    }//end of BeforePbc
    
    private static class AfterPbc implements IntegratorIntervalListener, java.io.Serializable {
        AfterPbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public int getPriority() {return 200;}//PBC is 100-199
        public void intervalAction(IntegratorIntervalEvent evt) {
            meter.iterator.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            //store last coordinate after pbc applied
           while(meter.iterator.hasNext()) {meter.rLast[i++].E(meter.iterator.nextAtom().coord.position());}
        }//end of intervalAction
        final MeterMeanSquareDisplacement meter;
    }//end of AfterPbc
    
    private int nAtoms = 0;
    AtomIteratorPhaseDependent iterator;
    Integrator integrator;
    private Vector[] rAccum, rLast;
    private final Space space;

}//end of class
