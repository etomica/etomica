/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.atom.IAtomList;
import etomica.atom.IAtom;
import etomica.config.IConformation;
import etomica.space.Space;

/*
 *  Geometry of Published Paracetamol Molecule (Orthorhombic)
 * 
 * @author Tai Tan
 */

public class ConformationParacetamolOrthorhombic implements IConformation {
	
	public ConformationParacetamolOrthorhombic(Space space) {
		this.space = space;
	}
	
	public void initializePositions(IAtomList list){
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		
		IAtom c1 = list.getAtom(SpeciesParacetamol.indexC[0]);
		x =   0.11531;
		y =   0.08253;
		z =   0.04163;
		c1.getPosition().E(new double [] {x, y, z});
		
		IAtom c2 = list.getAtom(SpeciesParacetamol.indexC[1]);
		x =   0.10650;
		y = - 1.18938;
		z = - 0.52643;
		c2.getPosition().E(new double [] {x, y, z});
		
		IAtom h1 = list.getAtom(SpeciesParacetamol.indexH[0]);
		x = - 0.79998;
		y = - 1.58600;
		z = - 0.93135;
		h1.getPosition().E(new double [] {x, y, z});
		
		IAtom c3 =  list.getAtom(SpeciesParacetamol.indexC[2]);
		x =   1.26324;
		y = - 1.94300;
		z = - 0.55501;
		c3.getPosition().E(new double [] {x, y, z});
		
		IAtom h2 = list.getAtom(SpeciesParacetamol.indexH[1]);
		x =   1.25933;
		y = - 2.92500;
		z = - 0.99131;
		h2.getPosition().E(new double [] {x, y, z});
		
		IAtom c4 = list.getAtom(SpeciesParacetamol.indexC[3]);
		x =   2.44760;
		y = - 1.44973;
		z = - 0.03026;
		c4.getPosition().E(new double [] {x, y, z});
		
		IAtom o1 = list.getAtom(SpeciesParacetamol.indexO[0]);
		x =   3.54566;
		y = - 2.23898;
		z = - 0.09715;
		o1.getPosition().E(new double [] {x, y, z});
		
		IAtom h5 = list.getAtom(SpeciesParacetamol.indexHp[0]);
		x =   4.29280;
		y = - 1.80970;
		z =   0.28439;
		h5.getPosition().E(new double [] {x, y, z});
		
		IAtom c5 = list.getAtom(SpeciesParacetamol.indexC[4]);
		x =   2.46224;
		y = - 0.18782;
		z =   0.53424;
		c5.getPosition().E(new double [] {x, y, z});
		
		IAtom h3 = list.getAtom(SpeciesParacetamol.indexH[2]);
		x =   3.37076;
		y =   0.21390;
		z =   0.95074;
		h3.getPosition().E(new double [] {x, y, z});
		
		IAtom c6 = list.getAtom(SpeciesParacetamol.indexC[5]);
		x =   1.30074;
		y =   0.56616;
		z =   0.57085;
		c6.getPosition().E(new double [] {x, y, z});
		
		IAtom h4 = list.getAtom(SpeciesParacetamol.indexH[3]);
		x =   1.32879;
		y =   1.54416;
		z =   1.02133;
		h4.getPosition().E(new double [] {x, y, z});
		
		IAtom n1 = list.getAtom(SpeciesParacetamol.indexN[0]);
		x = - 1.01685;
		y =   0.92837;
		z =   0.08726;
		n1.getPosition().E(new double [] {x, y, z});
		
		IAtom h6 = list.getAtom(SpeciesParacetamol.indexHp[1]);
		x = - 0.82779;
		y =   1.87097;
		z =   0.33401;
		h6.getPosition().E(new double [] {x, y, z});
		
		IAtom c7 = list.getAtom(SpeciesParacetamol.indexC[6]);
		x = - 2.32100;
		y =   0.61215;
		z = - 0.13549;
		c7.getPosition().E(new double [] {x, y, z});
		
		IAtom o2 = list.getAtom(SpeciesParacetamol.indexO[1]);
		x = - 2.71054;
		y = - 0.48505;
		z = - 0.41957;
		o2.getPosition().E(new double [] {x, y, z});
		
		IAtom c8 = list.getAtom(SpeciesParacetamol.indexC[7]);
		x = - 3.27190;
		y =   1.78724;
		z = - 0.03877;
		c8.getPosition().E(new double [] {x, y, z});
		
		IAtom h7 = list.getAtom(SpeciesParacetamol.indexH[4]);
		x = - 2.95649;
		y =   2.52790;
		z =   0.68746;
		h7.getPosition().E(new double [] {x, y, z});
		
		IAtom h8 = list.getAtom(SpeciesParacetamol.indexH[5]);
		x = - 4.25132;
		y =   1.41465;
		z =   0.22431;
		h8.getPosition().E(new double [] {x, y, z});
	
		IAtom h9 = list.getAtom(SpeciesParacetamol.indexC[6]);
		x = - 3.33710;
		y =   2.26662;
		z = - 1.01090;
		h9.getPosition().E(new double [] {x, y, z});
		
	}

	private static final long serialVersionUID = 1L;
    protected final Space space;
}
