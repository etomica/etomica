/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomType;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.Species;

/**
 *  
 * Species Naphthalene molecule
 * this is for TraPPE, the Naphthalene is rigid , LJ potential
 * reference: TraPPE 4, UA description of linear and branched alkanes and alkylbenzenes, Siepmann
 * 
 * @author shu
 * Oct, 20, 2010
 *
 */
public class SpeciesTraPPENaphthalene extends Species {

    public SpeciesTraPPENaphthalene(Space space) {
        this(space, false);
    }
    
    public SpeciesTraPPENaphthalene(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;
        
        chType = new AtomTypeLeaf(new ElementSimple("CH", 13.0107));
        cType = new AtomTypeLeaf(Carbon.INSTANCE);
        ////should change because it is not united atom!!!
        addChildType(chType);
        addChildType(cType);

        setConformation(new ConformationNaphthaleneTraPPE(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule Naphthalene = new Molecule(this, 10);
         // 2 Carbon without H, 8 Carbon with H
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Naphthalene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         conformation.initializePositions(Naphthalene.getChildList());
         return Naphthalene;
     }

     public IAtomType getCType() {
         return cType;
     }

     public IAtomType getCHType() {
         return chType;
     }


     public int getNumLeafAtoms() {
         return 10;
     }
    
    public final static int indexC1 = 0;
    public final static int indexC2 = 1;
    public final static int indexCH1 = 2;
    public final static int indexCH2 = 3;
    public final static int indexCH3 = 4;
    public final static int indexCH4 = 5;
    public final static int indexCH5 = 6;
    public final static int indexCH6 = 7;
    public final static int indexCH7 = 8;
    public final static int indexCH8 = 9;
       
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomTypeLeaf cType, chType;
}
