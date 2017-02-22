/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.config.IConformationOriented;
import etomica.space.IOrientation;
import etomica.space.ISpace;

 /**
  *  Conformation for HS Dimer
  *
  * @author Tai Boon Tan
  *
  */
public class ConformationHSDimer implements IConformationOriented, java.io.Serializable{
	
	public ConformationHSDimer(ISpace space, double L){
		this.space = space;
		this.L = L;
        orientationX = space.makeOrientation();
        IVectorMutable vectorX = space.makeVector();
        vectorX.setX(0, 1);
        orientationX.setDirection(vectorX);
	}

    public void initializePositions(IAtomList atomList, IOrientation orientation) {
        IVector orientationDir = orientation.getDirection();
        IAtom n1 = atomList.getAtom(SpeciesHSDimer.indexAtom1);
        n1.getPosition().Ea1Tv1(-0.5*L, orientationDir);
        
        IAtom n2 = atomList.getAtom(SpeciesHSDimer.indexAtom2);
        n2.getPosition().Ea1Tv1(+0.5*L, orientationDir);
    }

    public void initializePositions(IAtomList atomList) {
        initializePositions(atomList, orientationX);
	}
		
	protected final ISpace space;
	protected final double L;
	
	protected IOrientation orientationX;
	
	private static final long serialVersionUID = 1L;
}
