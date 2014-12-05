/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.Serializable;

import etomica.api.IAtom;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.space.ISpace;

/**
 * CoordinateDefinition implementation for monatomic molecules that are simply
 * leaf atoms. The class simply takes the u values to be real space
 * displacements from the nominal positions.
 * 
 * @author Andrew Schultz
 */

//when dealing with heterogeneous molecular systems we may need to introduce the mass in the definition

public class CoordinateDefinitionLeaf extends CoordinateDefinition implements
        Serializable {

    public CoordinateDefinitionLeaf(IBox box, Primitive primitive, ISpace space) {
        this(box, primitive, new BasisMonatomic(space), space);
    }
    
    public CoordinateDefinitionLeaf(IBox box, Primitive primitive, Basis basis, ISpace space) {
        super(box, space.D()*basis.getScaledCoordinates().length, primitive, basis, space);
        workVector = space.makeVector();
        u = new double[coordinateDim];
    }

    /**
     * Assigns the given array u to be the current position of the atom minus its lattice position
     */
    public double[] calcU(IMoleculeList atoms) {
        int j = 0;
        for (int i=0; i<atoms.getMoleculeCount(); i++) {
            IAtom a = atoms.getMolecule(i).getChildList().getAtom(0);
            IVectorMutable pos = a.getPosition();
            IVectorMutable site = getLatticePosition(a);
            workVector.Ev1Mv2(pos, site);
            for (int k=0; k<workVector.getD(); k++) {
                u[j+k] = workVector.getX(k);
            }
            j += workVector.getD();
        }
        return u;
    }

    public void initNominalU(IMoleculeList molecules) {
        //nothing to do -- lattice site is all information needed for u
    }

    /**
     * Sets the position of the atom to be its lattice position plus the offset u
     */
    public void setToU(IMoleculeList atoms, double[] newU) {
        int j = 0;
        for (int i=0; i<atoms.getMoleculeCount(); i++) {
            IAtom a = atoms.getMolecule(i).getChildList().getAtom(0);
            IVectorMutable pos = a.getPosition();
            for (int k=0; k<workVector.getD(); k++) {
                pos.setX(k, newU[j+k]);
            }
            j += workVector.getD();
            pos.PE(getLatticePosition(a));
        }
    }

    protected final IVectorMutable workVector;
    protected final double[] u;
    private static final long serialVersionUID = 1L;
}
