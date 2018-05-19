/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;

/**
 * Places atoms in a straight line.  Does not zero total momentum.
 * Can be used to implement a 1D linear conformation.
 *
 * @author David Kofke
 */

public class ConformationLinear implements IConformation, java.io.Serializable {
    
    public ConformationLinear(Space _space) {
        this(_space, 0.55);
    }
    public ConformationLinear(Space _space, double bondLength) {
    	this(_space, bondLength, makeDefaultAngles(_space));
    }
    
    private static double[] makeDefaultAngles(Space _space) {
        switch (_space.D()) {
            case 1: return new double[0];
            case 2: return new double[]{etomica.units.Degree.UNIT.toSim(45)};
            case 3: return new double[] {etomica.units.Degree.UNIT.toSim(45.), 0.0};
            default: throw new RuntimeException(_space.D()+" dimensional space?  I'm impressed.");
        }
    }
    
    public ConformationLinear(Space space, double bondLength, double[] initAngles) {
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
                setOrientation(Vector.of(new double[]{Math.cos(angle[0]), Math.sin(angle[0])}));
                return;
            case 3:
            	double[] ang = { Math.sin(angle[1])*Math.cos(angle[0]),
   			                     Math.sin(angle[1])*Math.sin(angle[0]),
   			                     Math.cos(angle[1]) };
                setOrientation(Vector.of(ang));
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

    public void initializePositions(IAtomList atomList) {
        int size = atomList.size();
        if(size == 0) return;

        double xNext = -bondLength*0.5*(size-1);
        int nLeaf = atomList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = atomList.get(iLeaf);
            a.getPosition().Ea1Tv1(xNext, orientation);
            xNext += bondLength;
        }
    }

    private static final long serialVersionUID = 1L;
    protected final Space space;
    protected double bondLength;
    private Vector orientation;
    private double[] angle;
}
