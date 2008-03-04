package etomica.config;
import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomActionTranslateTo;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * Places atoms in a straight line.
 * Can be used to implement a 1D linear conformation.
 *
 * @author David Kofke
 */

public class ConformationChainLinear extends ConformationChain {
    
    public ConformationChainLinear(ISimulation sim) {
        this(sim.getSpace(), 0.55);
    }
    public ConformationChainLinear(Space space, double bondLength) {
    	this(space, bondLength, new double[] {etomica.units.Degree.UNIT.toSim(45.), 0.0});
    }
    public ConformationChainLinear(Space space, double bondLength, double[] initAngles) {
        super(space);
        this.bondLength = bondLength;
        orientation = space.makeVector();
        angle = new double[space.D()];
        for(int i=0; i<initAngles.length; i++) setAngle(i,initAngles[i]);
    }
    public void setBondLength(double b) {
        bondLength = b;
    }

    public double getBondLength() {return bondLength;}
    public Dimension getBondLengthDimension() {return Length.DIMENSION;}
    
    //need to re-express this in terms of a Space.Orientation object
    public void setAngle(int i, double t) {//t in radians
        angle[i] = t;
        switch(angle.length) {
            case 1:
                return;
            case 2:
                setOrientation(new etomica.space2d.Vector2D(Math.cos(angle[0]),Math.sin(angle[0])));
                return;
            case 3:
                setOrientation(new etomica.space3d.Vector3D(Math.sin(angle[1])*Math.cos(angle[0]),
                                                  Math.sin(angle[1])*Math.sin(angle[0]),
                                                  Math.cos(angle[1])));
                return;
        }
    }
    public double getAngle(int i) {return angle[i];}
    public void setOrientation(IVector e) {orientation.E(e);}
    
    public void setOffset(IVector v) {
        orientation.E(v);
        bondLength = Math.sqrt(v.squared());
        orientation.TE(1.0/bondLength);
    }

    public void initializePositions(IAtomSet atomList) {
        int size = atomList.getAtomCount();
        if(size == 0) return;

        double xNext = -bondLength*0.5*(size-1);
        int nLeaf = atomList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomPositioned a = (IAtomPositioned)atomList.getAtom(iLeaf);
            a.getPosition().Ea1Tv1(xNext, orientation);
            xNext += bondLength;
        }
    }

    private static final long serialVersionUID = 1L;
    private double bondLength;
    private IVector orientation;
    private double[] angle;
	/* (non-Javadoc)
	 * @see etomica.ConformationChain#reset()
	 */
	protected void reset() {
		// TODO what happens here?
		
	}
	/* (non-Javadoc)
	 * @see etomica.ConformationChain#nextVector()
	 */
	protected IVector nextVector() {
		// TODO what happens here? 
		return null;
	}
}//end of ConformationLinear
      
