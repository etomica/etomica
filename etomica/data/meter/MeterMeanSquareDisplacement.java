package etomica.data.meter;
import etomica.EtomicaInfo;
import etomica.action.Action;
import etomica.api.IVector;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorBox;
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

public class MeterMeanSquareDisplacement extends DataSourceScalar {

    
    public MeterMeanSquareDisplacement(Space space, IntegratorBox integrator) {
        this(space, integrator, new AtomIteratorLeafAtoms());
    }

    public MeterMeanSquareDisplacement(Space space, IntegratorBox integrator, AtomIteratorBoxDependent iter) {
        super("Mean square displacement", Undefined.DIMENSION);
        this.space = space;
        this.integrator = integrator;
        setIterator(iter);
        BeforePbc beforePbc = new BeforePbc(this);
        integrator.addIntervalAction(beforePbc);
        integrator.setActionInterval(beforePbc, 50);
        AfterPbc afterPbc = new AfterPbc(this);
        integrator.addIntervalAction(afterPbc);
        integrator.setActionInterval(afterPbc, 200);
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Mean squared displacement, suitable for diffusion calculations");
        return info;
    }
    
    public void setIterator(AtomIteratorBoxDependent iterator) {
        if(iterator == null) {
            throw new NullPointerException("Cannot give a null iterator");
        }
        this.iterator = iterator;
        if(integrator.getBox() != null) {
            iterator.setBox(integrator.getBox());
            reset();
        } else {
            //throw an exception, because meter won't be informed when integrator has box set
            throw new IllegalStateException("Must first define Box for Integrator before constructing MeterMeanSquareDisplacement");
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
        for (IAtomPositioned a = (IAtomPositioned)iterator.nextAtom(); a != null;
             a = (IAtomPositioned)iterator.nextAtom()) {
            rAccum[i] = space.makeVector();
            rLast[i] = space.makeVector();
            rLast[i].E(a.getPosition());
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
    
    private static class BeforePbc implements Action, java.io.Serializable {
        BeforePbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public void actionPerformed() {
//            meter.iterator.setBox(meter.integrator.getBox()[0]);
            AtomIterator it = meter.iterator;
            it.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            for (IAtomPositioned a = (IAtomPositioned)it.nextAtom(); a != null;
                 a = (IAtomPositioned)it.nextAtom()) {
                IVector r = a.getPosition();
                meter.rAccum[i].PE(r);
                meter.rAccum[i].ME(meter.rLast[i]);
                meter.rLast[i].E(r);
                i++;
            }
        }
        private static final long serialVersionUID = 1L;
        final MeterMeanSquareDisplacement meter;
    }
    
    private static class AfterPbc implements Action, java.io.Serializable {
        AfterPbc(MeterMeanSquareDisplacement meter) {
            this.meter = meter;
        }
        public void actionPerformed() {
            AtomIterator it = meter.iterator;
            it.reset();
            int i = 0;
            //accumulate difference from last coordinate before pbc applied
            //store last coordinate after pbc applied
            for (IAtomPositioned a = (IAtomPositioned)it.nextAtom(); a != null;
                 a = (IAtomPositioned)it.nextAtom()) {
                meter.rLast[i++].E(a.getPosition());
            }
        }
        private static final long serialVersionUID = 1L;
        final MeterMeanSquareDisplacement meter;
    }
    
    private static final long serialVersionUID = 1L;
    private int nAtoms = 0;
    AtomIteratorBoxDependent iterator;
    IntegratorBox integrator;
    protected IVector[] rAccum, rLast;
    private final Space space;

}//end of class
