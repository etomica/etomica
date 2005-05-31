package etomica.space;

import etomica.Default;
import etomica.Space;
import etomica.lattice.IndexIteratorSequential;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polytope;
import etomica.math.geometry.Rectangle;
import etomica.math.geometry.Rectangular;

/**
 * Boundary that is in the shape of a rectangular parallelepiped.  
 * Periodicity in each direction is specified by subclass.
 */
/*
 * History Created on Jan 24, 2005 by kofke
 */
public abstract class BoundaryRectangular extends Boundary implements BoundaryPeriodic {

    public BoundaryRectangular(Space space, boolean[] periodicity) {
        super(space, makeShape(space));
        isPeriodic = (boolean[])periodicity.clone();
        dimensions = space.makeVector();
        dimensions.E(Default.BOX_SIZE);
        temp = space.makeVector();
        modShift = space.makeVector();
        dimensionsCopy = space.makeVector();
        dimensionsHalf = space.makeVector();
        indexIterator = new IndexIteratorSequential(space.D());
        needShift = new boolean[space.D()];//used by getOverflowShifts
        updateDimensions();
    }
    
    private static Polytope makeShape(Space space) {
        switch(space.D()) {
            case 1: return new LineSegment(space);
            case 2: return new Rectangle(space);
            case 3: return new Cuboid(space);
            default: throw new IllegalArgumentException("BoundaryRectangular not appropriate to given space");
        }
    }

    public final etomica.space.Vector dimensions() {
        return dimensionsCopy;
    }

    public etomica.space.Vector randomPosition() {
        temp.setRandom(dimensions);
        return temp;
    }

    private final void updateDimensions() {
        dimensionsHalf.Ea1Tv1(0.5, dimensions);
        dimensionsCopy.E(dimensions);
        ((Rectangular)shape).setEdgeLengths(dimensions);
    }


    public void setDimensions(etomica.space.Vector v) {
        dimensions.E(v);
        updateDimensions();
    }

    public double volume() {
        return shape.getVolume();
    }

    public boolean[] getPeriodicity() {
        return isPeriodic;
    }

    public double[][] imageOrigins(int nShells) {
        Vector workVector = space.makeVector();
        int shellFormula = (2 * nShells) + 1;
        int nImages = space.powerD(shellFormula) - 1;
        double[][] origins = new double[nImages][space.D()];
        indexIterator.setSize(shellFormula);
        indexIterator.reset();
        int k = 0;
        while(indexIterator.hasNext()) {
            int[] index = indexIterator.next();
            workVector.E(index);
            workVector.PE(-(double)nShells);
            if(workVector.isZero()) continue;
            workVector.TE(dimensions);
            workVector.assignTo(origins[k++]);
        }
        return origins;
    }

    public float[][] getOverflowShifts(Vector rr, double distance) {
       int D = space.D();
       int numVectors = 1;
       for (int i=1; i<D; i++) {
          if ((rr.x(i) - distance < 0.0) || (rr.x(i) + distance > dimensions.x(i))) {
             //each previous vector will need an additional copy in this dimension 
             numVectors *= 2;
             //remember that
             needShift[i] = true;
          }
          else {
             needShift[i] = false;
          }
       }
       
       if(numVectors == 1) return shift0;
       
       float[][] shifts = new float[numVectors][D];
       double[] rrArray = rr.toArray();
       for (int i=0; i<D; i++) {
          shifts[0][i] = (float)rrArray[i];
       }
       int iVector = 1;

       for (int i=0; iVector<numVectors; i++) {
          if (!needShift[i]) {
             //no shift needed for this dimension
             continue;
          }
          double delta = -dimensions.x(i);
          if (rr.x(i) - distance < 0.0) {
             delta = -delta;
          }
          //copy all previous vectors and apply a shift of delta to the copies
          for (int jVector=0; jVector<iVector; jVector++) {
             for (int j=0; j<D; j++) {
                 shifts[jVector+iVector][j] = shifts[jVector][j];
             }
             shifts[jVector+iVector][i] += delta;
          }
          iVector *= 2;
       }
       return shifts;
    }
    
    private final Vector temp;
    protected final Vector modShift;
    protected final Vector dimensions;
    protected final Vector dimensionsCopy;
    protected final Vector dimensionsHalf;
    private final IndexIteratorSequential indexIterator;
    private final boolean[] needShift;
    protected final boolean[] isPeriodic;
    protected final float[][] shift0 = new float[0][0];
    protected float[][] shift;

}
