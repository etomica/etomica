package etomica.config;
import etomica.action.AtomActionTranslateBy;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Dimension;
import etomica.units.Length;

/**
 * Places atoms in a straight line.  Does not zero total momentum.
 * Can be used to implement a 1D linear conformation.
 *
 * @author David Kofke
 */

public class ConformationLinear extends Conformation {
    
    public ConformationLinear(Simulation sim) {
        this(sim.space, 0.55*sim.getDefaults().atomSize);
    }
    public ConformationLinear(Space space, double bondLength) {
    	this(space, bondLength, new double[] {etomica.units.Degree.UNIT.toSim(45.), 0.0});
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
        atomIterator = new AtomIteratorArrayListSimple();
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
    public void setOrientation(Vector e) {orientation.E(e);}
    
    public void setOffset(Vector v) {
        orientation.E(v);
        bondLength = Math.sqrt(v.squared());
        orientation.TE(1.0/bondLength);
    }

    public void initializePositions(AtomArrayList atomList) {
        int size = atomList.size();
        if(size == 0) return;

        atomIterator.setList(atomList);
            
        double xNext = -bondLength*0.5*(size-1);
        atomIterator.reset();
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.nextAtom();
            if (!a.node.isLeaf()) {
                //initialize coordinates of child atoms
                Conformation config = a.type.creator().getConformation();
                config.initializePositions(((AtomTreeNodeGroup)a.node).childList);
            }
            moveToOrigin.actionPerformed(a);
            translationVector.Ea1Tv1(xNext, orientation);
            translator.actionPerformed(a);
            xNext += bondLength;
        }
    }

    private double bondLength;
    private Vector orientation;
    private double[] angle;
    private Vector translationVector;
    private AtomActionTranslateBy translator;
    private AtomActionTranslateTo moveToOrigin;
    private final AtomIteratorArrayListSimple atomIterator;
}//end of ConformationLinear
      
