/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.Molecule;
import etomica.chem.elements.Fluorine;
import etomica.chem.elements.Sulfur;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species SF6, 7 sites, LJ, rigid, no partial charge
 * Reference: Samios, Molecular force field investigation for sulfur hexafluoride: A computer simulation study
 * 
 * @author shu
 * 01-18-2013
 */
public class Species7SiteRigidSF6 extends Species {

    public final static int indexS = 0;
    public final static int indexF1 = 1;
    public final static int indexF2 = 2;
    public final static int indexF3 = 3;
    public final static int indexF4 = 4;
    public final static int indexF5 = 5;
    public final static int indexF6 = 6;
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final boolean isDynamic;
    protected final AtomType sType;
    protected final AtomType fType;
    public Species7SiteRigidSF6(Space space) {
        this(space, false);
    }
    public Species7SiteRigidSF6(Space space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;

        sType = new AtomType(Sulfur.INSTANCE);
        fType = new AtomType(Fluorine.INSTANCE);

        addChildType(sType);
        addChildType(fType);

        setConformation(new Conformation7SiteRigidSF6(space));
     }

     public IMolecule makeMolecule() {
         Molecule SF6 = new Molecule(this, 7);
         // 1 sulfur, 6 fluoride
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, sType) : new Atom(space, sType));
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, fType) : new Atom(space, fType));
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, fType) : new Atom(space, fType));
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, fType) : new Atom(space, fType));
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, fType) : new Atom(space, fType));
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, fType) : new Atom(space, fType));
         SF6.addChildAtom(isDynamic ? new AtomLeafDynamic(space, fType) : new Atom(space, fType));

         conformation.initializePositions(SF6.getChildList());
         return SF6;
     }

    public AtomType getSType() {
         return sType;
     }

    public AtomType getFType() {
         return fType;
     }

     public int getNumLeafAtoms() {
         return 7;
     }
     
    
}
