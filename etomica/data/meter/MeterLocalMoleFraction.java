package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorSpeciesDependent;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
import etomica.phase.Phase;
import etomica.species.Species;
import etomica.units.Dimension;

/**
 * Meter for measurement of species mole fraction within a specified subvolume
 */

public abstract class MeterLocalMoleFraction extends DataSourceScalar implements Meter
{
    public MeterLocalMoleFraction() {
        super("Local Mole Fraction",Dimension.FRACTION);
        setSpecies(null);
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
     * Accessor method to set which species mole-fraction or molar-density is averaged
     * To set to total number density, invoke with static ALL_SPECIES field as argument
     */
    public final void setSpecies(Species s) {species = s;}
    /**
     * Accessor method to get the current value of the species index
     *
     * @see #setSpeciesIndex
     */
    public final Species getSpecies() {return species;}
    
    /**
     * @return the current value of the local density or local mole fraction
     */
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        int totalSum = 0, speciesSum = 0;
        iterator.reset();
        while(iterator.hasNext()) {
            Atom m = iterator.nextAtom();
            if(shape.contains(m.coord.position())) {
                totalSum++;
                if(m.type.getSpecies() == species) speciesSum++;
            }
        }
        if(totalSum == 0) return Double.NaN;
        return (double)speciesSum/(double)totalSum;
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
            setShape(phase.boundary().getShape());
        }
    }

    /**
     * @return Returns the iterator.
     */
    public AtomIteratorSpeciesDependent getIterator() {
        return iterator;
    }
    /**
     * @param iterator The iterator to set.
     */
    public void setIterator(AtomIteratorSpeciesDependent iterator) {
        this.iterator = iterator;
    }

    private Phase phase;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private Species species;
    private AtomIteratorSpeciesDependent iterator = new AtomIteratorLeafAtoms();
    private Polytope shape;
}
