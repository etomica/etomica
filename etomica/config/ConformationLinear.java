package etomica.config;
import etomica.api.IAtomPositioned;
import etomica.api.IAtomSet;
import etomica.api.IConformation;
import etomica.api.IVector;
import etomica.space.ISpace;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * Places atoms in a straight line.  Does not zero total momentum.
 * Can be used to implement a 1D linear conformation.
 *
 * @author David Kofke
 */

public class ConformationLinear implements IConformation, java.io.Serializable {
    
    public ConformationLinear(ISpace _space) {
        this(_space, 0.55);
    }
    public ConformationLinear(ISpace _space, double bondLength) {
    	this(_space, bondLength, makeDefaultAngles(_space));
    }
    
    private static double[] makeDefaultAngles(ISpace _space) {
        switch (_space.D()) {
            case 1: return new double[0];
            case 2: return new double[]{etomica.units.Degree.UNIT.toSim(45)};
            case 3: return new double[] {etomica.units.Degree.UNIT.toSim(45.), 0.0};
            default: throw new RuntimeException(_space.D()+" dimensional space?  I'm impressed.");
        }
    }
    
    public ConformationLinear(ISpace space, double bondLength, double[] initAngles) {
        this.space = space;
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
    protected final ISpace space;
    protected double bondLength;
    private IVector orientation;
    private double[] angle;
}
