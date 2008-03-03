package etomica.atom;

import etomica.api.IVector;
import etomica.space.Space;


/**
 * Atom type for a rigid molecule with an orientation and (therefore) a moment
 * of inertia.  The molecule holds the orientation and the type holds the
 * moment.
 */
public class AtomTypeMoleculeOriented extends AtomTypeMolecule {

    public AtomTypeMoleculeOriented(Space space) {
        super(new AtomPositionCOM(space));
        moment = space.makeVector();
        this.space = space;
    }
    
    public void init() {
        // make a pretend molecule and calculate its moment of inertia
        IMolecule molecule = species.makeMolecule();
        AtomSet children = molecule.getChildList();
        conformation.initializePositions(children);
        IVector com = space.makeVector();
        com.E(positionDefinition.position(molecule));
        double[] I = new double[3];
        IVector xWork = space.makeVector();
        mass = 0;
        for (int i=0; i<children.getAtomCount(); i++) {
            IAtomPositioned atom = (IAtomPositioned)children.getAtom(i);
            xWork.Ev1Mv2(atom.getPosition(), com);
            double atomMass = ((AtomTypeLeaf)atom.getType()).getMass();
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

    /**
     * Returns the principle components of the moment of inertia of the
     * molecule within the body-fixed frame.  Do NOT modify the returned moment
     * of inertia returned.
     */
    public IVector getMomentOfInertia() {
        return moment;
    }
    
    public double getMass() {
        return mass;
    }

    private static final long serialVersionUID = 1L;
    protected final IVector moment;
    protected double mass;
    protected final Space space;
}