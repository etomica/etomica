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
  *  Reference paper: TraPPE: 4 UA description of linear and branched alkanes and alkylbenzenes, Siepmann et al
  *  modified from naphthalene
  *  kept the axis and coorinated unchanged
  * @author Shu Yang
  * 
  */
public class ConformationPhenanthreneTraPPE implements IConformation, java.io.Serializable{
	
	public ConformationPhenanthreneTraPPE(Space space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
		// put 2 C without H in y-axis , one is on the top, one is under the origin--this is naphthalene
		//Add 4 more cabon atoms
		
		IAtom n1 = atomList.get(SpeciesTraPPEPhenanthrene.indexC1);
		n1.getPosition().E(new double[] {0, a, 0});
		
		IAtom n2 = atomList.get(SpeciesTraPPEPhenanthrene.indexC2);
		n2.getPosition().E(new double[] {0, -a, 0});
						
		IAtom n3 = atomList.get(SpeciesTraPPEPhenanthrene.indexC3);
		n3.getPosition().E(new double[] {sqrt3a, -bond, 0});
				
		IAtom n4 = atomList.get(SpeciesTraPPEPhenanthrene.indexC4);
		n4.getPosition().E(new double[] {twosqrt3a, -a, 0});
		
		
		// put the other CH United atoms
		IAtom n5 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH1);
		n5.getPosition().E(new double[] {-sqrt3a, bond, 0});
		
		IAtom n6 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH2);
		n6.getPosition().E(new double[] {-twosqrt3a, a, 0});
		
		IAtom n7 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH3);
		n7.getPosition().E(new double[] {-twosqrt3a, -a, 0});
		
		IAtom n8 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH4);
		n8.getPosition().E(new double[] {-sqrt3a, -bond, 0});
		
		IAtom n9 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH5);
		n9.getPosition().E(new double[] {sqrt3a, -foura, 0});
		
		IAtom n10 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH6);
		n10.getPosition().E(new double[] {twosqrt3a, -fivea, 0});
		
		IAtom n11 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH7);
		n11.getPosition().E(new double[] {twosqrt3a, a, 0});
		
		IAtom n12 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH8);
		n12.getPosition().E(new double[] {sqrt3a, bond, 0});
		
		IAtom n13 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH9);
		
		n13.getPosition().E(new double[] {threesqrt3a, -foura, 0});
		
		IAtom n14 = atomList.get(SpeciesTraPPEPhenanthrene.indexCH10);
		n14.getPosition().E(new double[] {threesqrt3a, -bond, 0});
	
	}
	
    protected final Space space;
	protected static final double a = 0.7; // a is half of the bond length
	protected static final double bond = 1.4;
	protected static final double foura = a * 4;
	protected static final double fivea = a * 5;

	protected static final double sqrt3a = a * Math.sqrt(3);
	protected static final double twosqrt3a =  sqrt3a* 2;
	protected static final double threesqrt3a = sqrt3a* 3 ;
	protected static final double foursqrt3a = sqrt3a * 4 ;


	protected Vector vector;
	
	private static final long serialVersionUID = 1L;
}
