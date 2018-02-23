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
 *  Conformation of SF6
 *  7 LJ sites, rigid, no charge
 *  Reference: Samios, Molecular force field investigation for sulfur hexafluoride: A computer simulation study
 * 
 * @author shu
 * 01-18-2013
 */
public class Conformation7SiteRigidSF6 implements IConformation, java.io.Serializable{

	protected final Space space;
	protected final double bondL = 1.565;
	protected Vector vector;
	private static final long serialVersionUID = 1L;
	
	public Conformation7SiteRigidSF6(Space space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
		
		IAtom n1 = atomList.get(Species7SiteRigidSF6.indexS);
		n1.getPosition().E(new double[] {0, 0, 0});
		
		IAtom n2 = atomList.get(Species7SiteRigidSF6.indexF1);
		n2.getPosition().E(new double[] {0, bondL, 0});
		
		IAtom n3 = atomList.get(Species7SiteRigidSF6.indexF2);
		n3.getPosition().E(new double[] {0, -bondL, 0});
		
		IAtom n4 = atomList.get(Species7SiteRigidSF6.indexF3);
		n4.getPosition().E(new double[] {0, 0, bondL});
		
		IAtom n5 = atomList.get(Species7SiteRigidSF6.indexF4);
		n5.getPosition().E(new double[] {0, 0, -bondL});
		
		IAtom n6 = atomList.get(Species7SiteRigidSF6.indexF5);
		n6.getPosition().E(new double[] {bondL,0,0});
		
		IAtom n7 = atomList.get(Species7SiteRigidSF6.indexF6);
		n7.getPosition().E(new double[] {-bondL,0,0});
			
	}
	
 }
