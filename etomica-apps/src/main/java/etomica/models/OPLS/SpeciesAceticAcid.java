/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.OPLS;

import etomica.atom.IAtomType;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.space.Space;
import etomica.species.Species;

/**
 * Species for OPLS acetic acid
 * 
 * @author Hye Min Kim
 * Nov, 2011
 */
public class SpeciesAceticAcid extends Species {

    public SpeciesAceticAcid(Space space) {
    	
        super();
        
        this.space = space;
        
        cH3Type = new AtomTypeLeaf(new ElementSimple("cH3", 15.0107));//mass doesn't affect anything in MC simulation
        cType = new AtomTypeLeaf(Carbon.INSTANCE);
        dBOType = new AtomTypeLeaf(Oxygen.INSTANCE); 
        sBOType = new AtomTypeLeaf(Oxygen.INSTANCE);
        hType = new AtomTypeLeaf(Hydrogen.INSTANCE); 
        
        addChildType(cH3Type);
        addChildType(cType);
        addChildType(dBOType);
        addChildType(sBOType);
        addChildType(hType);

        setConformation(new ConformationAceticAcid(space)); 
     }

     public IMolecule makeMolecule() {
         Molecule aceticAcid = new Molecule(this, 3);

         aceticAcid.addChildAtom(new Atom(space, cH3Type));//0
         aceticAcid.addChildAtom(new Atom(space, cType));//1
         aceticAcid.addChildAtom(new Atom(space, dBOType));//2
         aceticAcid.addChildAtom(new Atom(space, sBOType));//3
         aceticAcid.addChildAtom(new Atom(space, hType));//4
 
         conformation.initializePositions(aceticAcid.getChildList());
         return aceticAcid;
     }
     
     public IAtomType getCH3Type() {
         return cH3Type;
     }
     public IAtomType getCType() {
         return cType;
     }
     
     public IAtomType getDBOType() {
         return dBOType;
     }
     
     public IAtomType getSBOType() {
         return sBOType;
     }

     public IAtomType getHType() {
         return hType;
     }

     public int getNumLeafAtoms() {
         return 5;
     }
    
    public final static int indexCH3  = 0;
    public final static int indexC  = 1;
    public final static int indexDBO   = 2;
    public final static int indexSBO   = 3;
    public final static int indexH  = 4;
  
    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected final AtomTypeLeaf cH3Type, cType, dBOType, sBOType, hType;
}
