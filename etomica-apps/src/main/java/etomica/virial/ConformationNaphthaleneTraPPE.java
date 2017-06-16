/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.config.IConformation;
import etomica.space.Space;

/**
  *  Conformation for Naphthalene
  *  Reference paper: TraPPE: 4 UA description of linear and branched alkanes and alkylbenzenes, Siepmann et al
  * 
 * @author Shu Yang
 *
 */
public class ConformationNaphthaleneTraPPE implements IConformation, java.io.Serializable{
	
	public ConformationNaphthaleneTraPPE(Space space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
			// put 2 C without H in y-axis , one is on the top, one is under the origin
		IAtom n1 = atomList.getAtom(SpeciesTraPPENaphthalene.indexC1);
		n1.getPosition().E(new double[] {0, halfofthebondlength, 0});
		
		IAtom n2 = atomList.getAtom(SpeciesTraPPENaphthalene.indexC2);
		n2.getPosition().E(new double[] {0, -halfofthebondlength, 0});
		
		// put the other CH United atoms
		IAtom n3 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH1);
		n3.getPosition().E(new double[] {-sqrt3ofhalfofthebondlength, bondlength, 0});
		
		IAtom n4 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH2);
		n4.getPosition().E(new double[] {-twicesqrt3ofhalfofthebondlength, halfofthebondlength, 0});
		
		IAtom n5 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH3);
		n5.getPosition().E(new double[] {-twicesqrt3ofhalfofthebondlength, -halfofthebondlength, 0});
		
		IAtom n6 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH4);
		n6.getPosition().E(new double[] {-sqrt3ofhalfofthebondlength, -bondlength, 0});
		
		IAtom n7 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH5);
		n7.getPosition().E(new double[] {sqrt3ofhalfofthebondlength, -bondlength, 0});
		
		IAtom n8 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH6);
		n8.getPosition().E(new double[] {twicesqrt3ofhalfofthebondlength, -halfofthebondlength, 0});
		
		IAtom n9 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH7);
		n9.getPosition().E(new double[] {twicesqrt3ofhalfofthebondlength, halfofthebondlength, 0});
		
		IAtom n10 = atomList.getAtom(SpeciesTraPPENaphthalene.indexCH8);
		n10.getPosition().E(new double[] {sqrt3ofhalfofthebondlength, bondlength, 0});
		
			
			
		
	}
	
    
	
	protected final Space space;
	protected static final double halfofthebondlength = 0.7;
	protected static final double bondlength = halfofthebondlength * 2;

	protected static final double sqrt3ofhalfofthebondlength = halfofthebondlength * Math.sqrt(3);
	protected static final double twicesqrt3ofhalfofthebondlength = sqrt3ofhalfofthebondlength * 2;
	protected Vector vector;
	
	private static final long serialVersionUID = 1L;
}
