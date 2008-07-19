package etomica.species;

import etomica.api.IAtomLeaf;
import etomica.api.IMolecule;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.atom.MoleculeOriented;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.chem.elements.ElementSimple;
import etomica.space.ISpace;

/**
 * Species in which molecules are made of a single atom.  The molecule itself
 * holds the orientation.
 *
 * @author Andrew Schultz
 */
public class SpeciesSpheresRotatingMolecule extends SpeciesSpheresMono implements ISpeciesOriented {
    
    public SpeciesSpheresRotatingMolecule(ISimulation sim, ISpace _space) {
        this(sim, _space, makeNominalMoment(_space));
    }

    protected static final IVector makeNominalMoment(ISpace space) {
        IVector m = space.makeVector();
        m.E(1);
        return m;
    }

    public SpeciesSpheresRotatingMolecule(ISimulation sim, ISpace _space, IVector moment) {
        super(_space, sim.isDynamic(), new AtomTypeSphere(new ElementSimple(sim), 1.0));
        this.moment = moment;
    }

    /**
     * Constructs a new group.
     */
     public IMolecule makeMolecule() {
         MoleculeOriented group = isDynamic ? new MoleculeOrientedDynamic(space, this)
                                            : new MoleculeOriented(space, this);
         group.addChildAtom(makeLeafAtom());
         return group;
     }

    protected IAtomLeaf makeLeafAtom() {
        return new AtomLeaf(space, leafAtomType);
    }

    public double getMass() {
        return leafAtomType.getMass();
    }

    public IVector getMomentOfInertia() {
        return moment;
    }

    protected IVector moment;
    private static final long serialVersionUID = 1L;
}
