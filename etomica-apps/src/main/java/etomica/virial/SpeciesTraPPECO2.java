/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.Molecule;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * 
 * 
 * Species CO2 molecule
 *this is for TraPPE, the CO2 is rigid , LJ potential + QQ 
 * 
 * @author Shu Yang
 * Oct, 20, 2010
 *
 */
public class SpeciesTraPPECO2 extends Species {

    public final static int indexC = 0;
    public final static int indexOleft = 1;
    public final static int indexOright = 2;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomType cType, oType;
    public SpeciesTraPPECO2(Space space) {
        this(space, false);
    }
    public SpeciesTraPPECO2(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;

        cType = new AtomType(Carbon.INSTANCE);
        oType = new AtomType(Oxygen.INSTANCE);
        addChildType(cType);
        addChildType(oType);

        setConformation(new ConformationCO2(space));
     }

    public IMolecule makeMolecule() {
         Molecule CO2 = new Molecule(this, 3);
         CO2.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         CO2.addChildAtom(isDynamic ? new AtomLeafDynamic(space, oType) : new Atom(space, oType));
         CO2.addChildAtom(isDynamic ? new AtomLeafDynamic(space, oType) : new Atom(space, oType));

         conformation.initializePositions(CO2.getChildList());
         return CO2;
     }

    public AtomType getCarbonType() {
         return cType;
     }

    public AtomType getOxygenType() {
         return oType;
     }

     public int getNumLeafAtoms() {
         return 3;
     }
}
