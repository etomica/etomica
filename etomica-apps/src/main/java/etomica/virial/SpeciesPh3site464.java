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
import etomica.chem.elements.ElementSimple;
import etomica.space.Space;
import etomica.species.Species;

/**
 *  
 * Species Phenanthrene molecule
 * rigid , LJ potential, no charge, 3 site model.each benzene ring is a site, 464 model
 * reference: Iwai, monte carlo sim of Naphthalene, phenathlene,anthracene in SCF 1998
 * @author shu
 * March, 9, 2011
 *
 */
public class SpeciesPh3site464 extends Species {

    public SpeciesPh3site464(Space space) {
        this(space, false);
    }
    
    public SpeciesPh3site464(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;
        // the virial simulation  doesnt care this , but still I added this info here FYI
        //CMass is the middle site, including 6 carbons and 6 hydrogens
        //CHMass is the side site, including 4 carbons and 4 hydrogens
        double CMass = 6 * 12.0107 + 1 * 2;
        double CHMass = 4 * 12.0107 + 4 * 5;
        chType = new AtomTypeLeaf(new ElementSimple("CH", CHMass));
        cType = new AtomTypeLeaf(new ElementSimple("C", CMass));
               
        //should change because it is not united atom!!!
        addChildType(chType);
        addChildType(cType);

        setConformation(new ConformationPh3site(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule Phenanthrene = new Molecule(this, 10);

         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, cType) : new Atom(space, cType));
         
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         Phenanthrene.addChildAtom(isDynamic ? new AtomLeafDynamic(space, chType) : new Atom(space, chType));
         conformation.initializePositions(Phenanthrene.getChildList());
         return Phenanthrene;
     }

     public IAtomType getCType() {
         return cType;
     }

     public IAtomType getCHType() {
         return chType;
     }

        // return the total number of the atoms
     public int getNumLeafAtoms() {
         return 3;
     }
    
    public final static int indexC1 = 0;
    public final static int indexCH1 = 1;
    public final static int indexCH2 = 2;
   
    
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomTypeLeaf cType, chType;
}
