/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;


import etomica.atom.Atom;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.chem.elements.ElementSimple;
import etomica.molecule.IMolecule;
import etomica.molecule.Molecule;
import etomica.space.Space;
import etomica.species.Species;
/**
 * Species for TraPPE-Explicit hydrogen,Siepmann
 * 
 * @shu
 * 01-30-2013
 */
  public class SpeciesAlkaneEH extends Species  {

	  	private static final long serialVersionUID = 1L;
	    protected final Space space;
	    protected final int numCarbons;
    protected final int totalAtoms;
    public int indexC_3_1;//CH3 on the left
	    public int indexC_3_2;//CH3 on the right
	    public int numH;
    //    protected final AtomType c_3Type, c_2Type, hType;
    public int numCH2;
	    protected boolean isDynamic;
	    public SpeciesAlkaneEH(Space space, int numCarbons) {
	    	super();
	        this.space = space;
	        this.numCarbons = numCarbons;
	        indexC_3_1 = 0;
	        indexC_3_2 = numCarbons - 1 ;
	        
	        numCH2 = numCarbons -2;
	        numH = numCarbons * 2 + 2;
	        totalAtoms = numCarbons * 3 + 2;

            addChildType(new AtomType(new ElementSimple("c_3", 1.0)));
            addChildType(new AtomType(new ElementSimple("c_2", 1.0)));
            addChildType(new AtomType(new ElementSimple("h", 1.0)));

            setConformation(new ConformationAlkaneEH(space, this));
	     }
	
	    public IMolecule makeMolecule() {
	    	
	    	Molecule molecule = new Molecule(this, totalAtoms);
	    	
	    	molecule.addChildAtom(makeLeafAtom(childTypes[0])); // the left C of CH3, global index [0]
	    	
	       	for(int j = 0; j < numCH2; j++) {//CH2, global index: [1]~[n-2]
	    		molecule.addChildAtom(makeLeafAtom(childTypes[1]));
	       	}
	    	molecule.addChildAtom(makeLeafAtom(childTypes[0])); // the right C of CH3,global index [n-1]
	
	    	for(int i = 0;  i <  numH;  i++) {//H 
	    		molecule.addChildAtom(makeLeafAtom(childTypes[2]));
	    	}
	    	conformation.initializePositions(molecule.getChildList());
			return molecule;
	    }

    public AtomType getC_3Type() {
        return childTypes[0];
	     }

    public AtomType getC_2Type() {
        return childTypes[1];
	     }

    public AtomType getHType() {
        return childTypes[2];
	     }
	     public int getNumLeafAtoms() {
	         return totalAtoms;
	     }

    protected IAtom makeLeafAtom(AtomType leafType) {
        return isDynamic ? new AtomLeafDynamic(space, leafType) : new Atom(space, leafType);
	     }
	
	}
