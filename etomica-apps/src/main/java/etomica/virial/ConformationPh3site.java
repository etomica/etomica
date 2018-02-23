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
  *  Conformation for Phenanthrene
  *  Reference paper: Iwai
  * 3 site group H to the two side benzene rings , eihter 464 or 545 model
  * modified from ConformationAnthracene
 * @author Shu
 * March,9,2011
 */
public class ConformationPh3site implements IConformation, java.io.Serializable{
	
	public ConformationPh3site(Space space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
			
		IAtom n1 = atomList.get(SpeciesPh3site545.indexC1);
		n1.getPosition().E(new double[] {0, 0, 0});
		
		IAtom n2 = atomList.get(SpeciesAnthracene3site545.indexCH1);
		n2.getPosition().E(new double[] {-bondlength, 0, 0});
		
		IAtom n3 = atomList.get(SpeciesAnthracene3site545.indexCH2);
		n3.getPosition().E(new double[] {half_bondlength, minus_half_sqrt_bondlength, 0});
				
	}
		
	protected final Space space;
	protected static final double bondlength = 2.42;//converst nm to angstrom
	protected static final double half_bondlength = bondlength / 2;
	protected static final double minus_half_sqrt_bondlength = - Math.sqrt(3) / 2 * bondlength ; 
	protected Vector vector;
	
	private static final long serialVersionUID = 1L;
}
