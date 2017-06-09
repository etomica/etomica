/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.atom.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.Molecule;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Nitrogen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * 
 * 
 * Species nitrogen molecule (shell model) 
 * 
 * Reference: Fabianski R. et al, Calculations on the stability of low temperature solid nitrogen
 *             phases, JCP 112(15) 6745 (2000)
 *             
 * @author Tai Boon Tan
 *
 */
public class SpeciesN2ShellModel extends Species {

    public final static int indexN1 = 0;
    public final static int indexN2 = 1;
    public final static int indexCenter = 2;
    public final static int indexP1left = 3;
    public final static int indexP1right = 4;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomType nType, pType;
    public SpeciesN2ShellModel(Space space) {
        this(space, false);
    }
    public SpeciesN2ShellModel(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;

        nType = new AtomType(Nitrogen.INSTANCE);
        pType = new AtomType(new ElementSimple("P", 1.0));
        addChildType(nType);
        addChildType(pType);

        setConformation(new ConformationNitrogenShellModel(space));
     }

     public IMolecule makeMolecule() {
         Molecule nitrogenShellModel = new Molecule(this, 5);
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, nType) : new Atom(space, nType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, nType) : new Atom(space, nType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));
         nitrogenShellModel.addChildAtom(isDynamic ? new AtomLeafDynamic(space, pType) : new Atom(space, pType));

         conformation.initializePositions(nitrogenShellModel.getChildList());
         return nitrogenShellModel;
     }

    public AtomType getNitrogenType() {
         return nType;
     }

    public AtomType getPType() {
         return pType;
     }

     public int getNumLeafAtoms() {
         return 5;
     }
}
