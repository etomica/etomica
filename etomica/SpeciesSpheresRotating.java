package etomica;

/**
 * Species in which molecules are made of a single atom of type OrientedSphere
 *
 * @author David Kofke
 * @see AtomType.OrientedSphere
 * 
 */
public class SpeciesSpheresRotating extends SpeciesDisks {
    
    public SpeciesSpheresRotating() {
        this(Simulation.instance);
    }
    public SpeciesSpheresRotating(Simulation sim) {
        this(sim, Default.MOLECULE_COUNT);
    }
    public SpeciesSpheresRotating(int nM) {
        this(Simulation.instance, nM);
    }
    public SpeciesSpheresRotating(Simulation sim, int nM) {
        super(sim, nM, 1, new AtomType.OrientedSphere(Default.ATOM_MASS, Default.ATOM_COLOR, Default.ATOM_SIZE));
    }
}


