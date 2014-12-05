/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

	package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.api.IVectorMutable;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialCalculation;
import etomica.space.Space;

/**
 * 
 * @author Tai Boon Tan
 */
public class PCEnergySumEinsteinCrystalExpansion implements PotentialCalculation, java.io.Serializable {

    /**
     * 
     * 
	 */
	public PCEnergySumEinsteinCrystalExpansion(){
		dr = Space.makeVector(3);
		latticeDistance = Space.makeVector(3);
	
	}
	
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
			
			atom0 = atoms.getAtom(0);
			atom1 = atoms.getAtom(1);
			
			dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
			box.getBoundary().nearestImage(dr);
			
			double distanceScalar = dr.squared();
			double rdu = ((Potential2SoftSpherical)potential).du(distanceScalar);
			   
			latticeDistance.Ev1Mv2(initialLatticePosition[atom1.getLeafIndex()], initialLatticePosition[atom0.getLeafIndex()]);
			box.getBoundary().nearestImage(latticeDistance);
			
			sum += rdu*(dr.dot(latticeDistance))/distanceScalar;
		
	}
	
	/**
	 * Sets the energy sum to zero, typically to begin a new energy-sum calculation.
	 */
	public void zeroSum() {
		sum = 0.0;
	}

	/**
	 * Returns the current value of the energy sum.
	 */
	public double getSum() {
        return sum;
    }
	
	public IVectorMutable[] getInitialLatticePosition() {
		return initialLatticePosition;
	}

	public void setInitialLatticePosition(IVectorMutable[] initialLatticePosition) {
		this.initialLatticePosition = initialLatticePosition;
	}

	public IBox getBox() {
		return box;
	}

	public void setBox(IBox box) {
		this.box = box;
	}
	
    private static final long serialVersionUID = 1L;
	protected  double sum = 0.0;
	public IAtom atom0, atom1;
	protected IVectorMutable[] initialLatticePosition;
	protected IVectorMutable latticeDistance, dr;
	protected IBox box;

}
