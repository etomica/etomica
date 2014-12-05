/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.IAtomType;
import etomica.api.IMolecule;
import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.Molecule;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Oxygen;
import etomica.models.nitrogen.ConformationNitrogen;
import etomica.space.ISpace;
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

    public SpeciesTraPPECO2(ISpace space) {
        this(space, false);
    }
    
    public SpeciesTraPPECO2(ISpace space, boolean isDynamic) {
        super();
        this.space = space;
        this.isDynamic = isDynamic;
        
        cType = new AtomTypeLeaf(Carbon.INSTANCE);
        oType = new AtomTypeLeaf(Oxygen.INSTANCE);
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

     public IAtomType getCarbonType() {
         return cType;
     }

     public IAtomType getOxygenType() {
         return oType;
     }


     public int getNumLeafAtoms() {
         return 3;
     }
    
    public final static int indexC = 0;
    public final static int indexOleft = 1;
    public final static int indexOright  = 2;
   
    
    private static final long serialVersionUID = 1L;
    protected final ISpace space;
    protected final boolean isDynamic;
    protected final AtomTypeLeaf cType, oType;
}