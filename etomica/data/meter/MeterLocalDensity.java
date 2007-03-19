package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
import etomica.phase.Phase;
import etomica.units.DimensionRatio;
import etomica.units.Quantity;
import etomica.units.Volume;

/**
 * Meter for measurement of density within a specified subvolume
 */
 
public abstract class MeterLocalDensity extends DataSourceScalar {
    public MeterLocalDensity() {
        super("Local Density",new DimensionRatio(Quantity.DIMENSION, Volume.DIMENSION));
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Local number density in a subregion of a phase");
        return info;
    }
    
    public void setShape(Polytope shape) {
        this.shape = shape;
    }

    public Polytope getShape() {
        return shape;
    }

    /**
     * @return the current value of the local density or local mole fraction
     */
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        //compute local molar density
        int nSum = 0;
        iterator.reset();
        while(iterator.hasNext()) {
            if(shape.contains(((AtomLeaf)iterator.nextAtom()).getCoord().getPosition())) nSum++;
        }
        return nSum/shape.getVolume();
    }
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        iterator.setPhase(phase);
        if (shape == null) {
            setShape(phase.getBoundary().getShape());
        }
    }

    /**
     * @return Returns the iterator.
     */
    public AtomIteratorPhaseDependent getIterator() {
        return iterator;
    }
    /**
     * @param iterator The iterator to set.
     */
    public void setIterator(AtomIteratorPhaseDependent iterator) {
        this.iterator = iterator;
    }

    private Phase phase;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private AtomIteratorPhaseDependent iterator = new AtomIteratorLeafAtoms();
    private Polytope shape;
}
