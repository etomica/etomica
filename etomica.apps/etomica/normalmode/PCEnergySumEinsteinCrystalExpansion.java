	package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.api.IPotential;
import etomica.api.IVector;
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
	
	public void doCalculation(IAtomSet atoms, IPotential potential) {
			
			atom0 = atoms.getAtom(0);
			atom1 = atoms.getAtom(1);
			
			dr.Ev1Mv2(((IAtomPositioned)atom1).getPosition(), ((IAtomPositioned)atom0).getPosition());
			box.getBoundary().nearestImage(dr);
			
			double distanceScalar = dr.squared();
			double rdu = ((Potential2SoftSpherical)potential).du(distanceScalar);
			   
			latticeDistance.Ev1Mv2(initialLatticePosition[((IAtomLeaf)atom1).getLeafIndex()], initialLatticePosition[((IAtomLeaf)atom0).getLeafIndex()]);
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
	
	public IVector[] getInitialLatticePosition() {
		return initialLatticePosition;
	}

	public void setInitialLatticePosition(IVector[] initialLatticePosition) {
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
	protected IVector[] initialLatticePosition;
	protected IVector latticeDistance, dr;
	protected IBox box;

}
