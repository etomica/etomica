package etomica.species;

import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IMolecule;
import etomica.api.IVector;
import etomica.atom.AtomPositionCOM;
import etomica.space.Space;


/**
 * Atom type for a rigid molecule with an orientation and (therefore) a moment
 * of inertia.  The molecule holds the orientation and the type holds the
 * moment.
 */
public abstract class SpeciesOriented extends Species implements ISpeciesOriented {

    public SpeciesOriented(Space space) {
        super(new AtomPositionCOM(space));
        moment = space.makeVector();
        this.space = space;
        // we could call init here, but the subclass might not be ready to make
        // a molecule yet.  Subclass will call init when it's ready.
    }
    
    protected void init() {
        // make a pretend molecule and calculate its moment of inertia
        IMolecule molecule = makeMolecule();
        IAtomSet children = molecule.getChildList();
        conformation.initializePositions(children);
        IVector com = space.makeVector();
        com.E(positionDefinition.position(molecule));
        double[] I = new double[3];
        IVector xWork = space.makeVector();
        mass = 0;
        for (int i=0; i<children.getAtomCount(); i++) {
            IAtomPositioned atom = (IAtomPositioned)children.getAtom(i);
            xWork.Ev1Mv2(atom.getPosition(), com);
            double atomMass = ((IAtomTypeLeaf)atom.getType()).getMass();
            mass += atomMass;
            for (int j=0; j<3; j++) {
                for (int k=0; k<3; k++) {
                    if (j==k) continue;
                    I[j] += atomMass*xWork.x(k)*xWork.x(k);
                }
            }
        }
        moment.E(I);
    }

    /* (non-Javadoc)
     * @see etomica.atom.ISpeciesOriented#getMomentOfInertia()
     */
    public IVector getMomentOfInertia() {
        return moment;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.ISpeciesOriented#getMass()
     */
    public double getMass() {
        return mass;
    }

    private static final long serialVersionUID = 1L;
    protected final IVector moment;
    protected double mass;
    protected final Space space;
}