package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceScalar;
import etomica.math.geometry.Polytope;
import etomica.phase.Phase;
import etomica.space.IVector;
import etomica.species.Species;
import etomica.units.Fraction;

/**
 * Meter for measurement of species mole fraction within a specified subvolume
 */
public class MeterLocalMoleFraction extends DataSourceScalar {

    public MeterLocalMoleFraction() {
        super("Local Mole Fraction",Fraction.DIMENSION);
        setSpecies(null);
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Local number density in a subregion of a phase");
        return info;
    }

    /**
     * Sets the subvolume shape for the mole fraction calculation.
     */
    public void setShape(Polytope shape) {
        this.shape = shape;
    }

    /**
     * Returns the subvolume shape for the mole fraction calculation.
     */
    public Polytope getShape() {
        return shape;
    }
    
    /**
     * Sets the origin of the subvolume (Polytopes typically have their center
     * at 0).  The shape origin is 0 by default.
     */
    public void setShapeOrigin(IVector newShapeOrigin) {
        shapeOrigin = newShapeOrigin;
    }
    
    /**
     * Returns the origin of the subvolume.
     */
    public IVector getShapeOrigin() {
        return shapeOrigin;
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
            AtomLeaf a = (AtomLeaf)iterator.nextAtom();
            tempVec.Ev1Mv2(a.getPosition(), shapeOrigin);
            if(shape.contains(tempVec)) {
                totalSum++;
                if(a.getType().getSpecies() == species) speciesSum++;
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
    public void setPhase(Phase newPhase) {
        phase = newPhase;
        tempVec = phase.getSpace().makeVector();
        shapeOrigin = phase.getSpace().makeVector();
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

    private static final long serialVersionUID = 1L;
    private Phase phase;
    /**
     * Class variable used to specify that all species are included in number-density calculation
     */
    private Species species;
    private AtomIteratorPhaseDependent iterator = new AtomIteratorLeafAtoms();
    private Polytope shape;
    private IVector shapeOrigin;
    private IVector tempVec;
}
