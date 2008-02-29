package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.api.IBox;
import etomica.atom.IAtomPositioned;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorBoxDependent;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
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
        EtomicaInfo info = new EtomicaInfo("Local number density in a subregion of a box");
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
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        //compute local molar density
        int nSum = 0;
        iterator.reset();
        for (IAtomPositioned atom = (IAtomPositioned)iterator.nextAtom(); atom != null;
             atom = (IAtomPositioned)iterator.nextAtom()) {
            if(shape.contains(atom.getPosition())) nSum++;
        }
        return nSum/shape.getVolume();
    }
    
    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
        iterator.setBox(box);
        if (shape == null) {
            setShape(box.getBoundary().getShape());
        }
    }

    /**
     * @return Returns the iterator.
     */
    public AtomIteratorBoxDependent getIterator() {
        return iterator;
    }
    /**
     * @param iterator The iterator to set.
     */
    public void setIterator(AtomIteratorBoxDependent iterator) {
        this.iterator = iterator;
    }

    private IBox box;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private AtomIteratorBoxDependent iterator = new AtomIteratorLeafAtoms();
    private Polytope shape;
}
