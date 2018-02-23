/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.species;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePositionCOM;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Atom type for a rigid molecule with an orientation and (therefore) a moment
 * of inertia.  The molecule holds the orientation and the type holds the
 * moment.
 */
public abstract class SpeciesOriented extends Species implements ISpeciesOriented {

    public SpeciesOriented(Space space) {
        super();
        moment = space.makeVector();
        this.space = space;
        // we could call init here, but the subclass might not be ready to make
        // a molecule yet.  Subclass will call init when it's ready.
    }
    
    protected void init() {
        // make a pretend molecule and calculate its moment of inertia
        IMolecule molecule = makeMolecule();
        IAtomList children = molecule.getChildList();
        conformation.initializePositions(children);
        Vector com = space.makeVector();
        MoleculePositionCOM positionCOM = new MoleculePositionCOM(space);
        com.E(positionCOM.position(molecule));
        double[] I = new double[3];
        Vector xWork = space.makeVector();
        mass = 0;
        for (int i = 0; i<children.size(); i++) {
            IAtom atom = children.get(i);
            xWork.Ev1Mv2(atom.getPosition(), com);
            double atomMass = atom.getType().getMass();
            mass += atomMass;
            for (int j=0; j<3; j++) {
                for (int k=0; k<3; k++) {
                    if (j==k) continue;
                    I[j] += atomMass*xWork.getX(k)*xWork.getX(k);
                }
            }
        }
        moment.E(I);
    }

    /* (non-Javadoc)
     * @see etomica.atom.ISpeciesOriented#getMomentOfInertia()
     */
    public Vector getMomentOfInertia() {
        return moment;
    }
    
    /* (non-Javadoc)
     * @see etomica.atom.ISpeciesOriented#getMass()
     */
    public double getMass() {
        return mass;
    }

    private static final long serialVersionUID = 1L;
    protected final Vector moment;
    protected double mass;
    protected final Space space;
}
