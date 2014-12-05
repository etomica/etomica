/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.api.IAtom;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class BiasVolumeCube extends BiasVolume {
    
    private static final long serialVersionUID = 1L;
    private final IVectorMutable dimensions;
    private final IVectorRandom work;
    private final IRandom random;
    private IBoundary boundary;
    
    public BiasVolumeCube(ISpace space, IRandom random){
        super(space);
        this.random = random;
        dimensions = space.makeVector();
        work = (IVectorRandom)space.makeVector();
        dimensions.E(2.0);//size of the cube
    }
    
    public void setBiasCubeDimensions(IVector dim) {
        dimensions.E(dim);
    }
    public void setBiasCubeDimensions(double d) {dimensions.E(d);}
    public IVector getBiasCubeDimensions() {return dimensions;}
    
    public void setBox(IBox box) {
    	boundary = box.getBoundary();
    }
    
    public double biasVolume() {
        double prod = 1.0;
        for(int i=0; i<space.D(); i++) {
            prod *= dimensions.getX(i);
        }
        return prod;
    }
    
    
    // Insert atom1 in to the Bonding region of atom2 
    //Bonding region is a cube /rectangle for 2 D
    public void biasInsert(IAtom atom1, IAtom atom2) {
        work.setRandomCube(random);
        work.TE(dimensions);
        atom1.getPosition().Ev1Pv2(atom2.getPosition(), work);
    }

    /**
     *Function to check for bonding
     */
    
    public boolean isAssociated(IAtom atom1, IAtom atom2){
    
        work.E(atom2.getPosition());
        work.ME(atom1.getPosition());
        boundary.nearestImage(work);
        work.DE(dimensions);
        for (int i = 0; i < work.getD(); i++) {
        	if (Math.abs(work.getX(i)) > 0.5 ) {
        		return false;
        	}
        }
        return true;
    }
}
