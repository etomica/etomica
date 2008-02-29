package etomica.config;
import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomActionTranslateTo;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.simulation.ISimulation;
import etomica.space.Space;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * Places atoms in a straight line.  Does not zero total momentum.
 * Can be used to implement a 1D linear conformation.
 *
 * @author David Kofke
 */

public class ConformationLinear extends Conformation {
    
    public ConformationLinear(ISimulation sim) {
        this(sim.getSpace(), 0.55);
    }
    public ConformationLinear(Space space, double bondLength) {
    	this(space, bondLength, makeDefaultAngles(space));
    }
    
    private static double[] makeDefaultAngles(Space space) {
        switch (space.D()) {
            case 1: return new double[0];
            case 2: return new double[]{etomica.units.Degree.UNIT.toSim(45)};
            case 3: return new double[] {etomica.units.Degree.UNIT.toSim(45.), 0.0};
            default: throw new RuntimeException(space.D()+" dimensional space?  I'm impressed.");
        }
    }
    
    public ConformationLinear(Space space, double bondLength, double[] initAngles) {
        super(space);
        this.bondLength = bondLength;
        orientation = space.makeVector();
        angle = new double[space.D()];
        for(int i=0; i<initAngles.length; i++) setAngle(i,initAngles[i]);
        translator = new AtomActionTranslateBy(space);
        moveToOrigin = new AtomActionTranslateTo(space);
        translationVector = translator.getTranslationVector();
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

    public void initializePositions(AtomSet atomList) {
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
    protected double bondLength;
    private IVector orientation;
    private double[] angle;
    private IVector translationVector;
    private AtomActionTranslateBy translator;
    private AtomActionTranslateTo moveToOrigin;
}
