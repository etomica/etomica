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
import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.Species;

/**
 *  
 * Species Phenanthrene molecule
 * this is for TraPPE, the Phenanthrene is rigid , LJ potential
 * reference: TraPPE 4, UA description of linear and branched alkanes and alkylbenzenes, Siepmann
 * modified from species of phenanthrene C14H10
 * @author shu
 * March.19.2011
 *
 */
public class SpeciesTraPPEPhenanthrene extends Species {

    public final static int indexC1 = 0;
    public final static int indexC2 = 1;
    public final static int indexC3 = 2;
    public final static int indexC4 = 3;
    public final static int indexCH1 = 4;
    public final static int indexCH2 = 5;
    public final static int indexCH3 = 6;
    public final static int indexCH4 = 7;
    public final static int indexCH5 = 8;
    public final static int indexCH6 = 9;
    public final static int indexCH7 = 10;
    public final static int indexCH8 = 11;
    public final static int indexCH9 = 12;
    public final static int indexCH10 = 13;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomType cType, chType;
    public SpeciesTraPPEPhenanthrene(Space space) {
        this(space, false);
    }
    public SpeciesTraPPEPhenanthrene(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;

        chType = new AtomType(new ElementSimple("CH", 13.0107));
        cType = new AtomType(Carbon.INSTANCE);
        ////should change because it is not united atom!!!
        addChildType(chType);
        addChildType(cType);

        setConformation(new ConformationPhenanthreneTraPPE(space));
     }

    public IMolecule makeMolecule() {
         Molecule Phenanthrene = new Molecule(this, 14);
         // 4 Carbon without H, 8 Carbon with H
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));

         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));

         conformation.initializePositions(Phenanthrene.getChildList());
         return Phenanthrene;
     }

    public AtomType getCType() {
         return cType;
     }

    public AtomType getCHType() {
         return chType;
     }

     public int getNumLeafAtoms() {
         return 14;
     }
}
