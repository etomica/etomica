package etomica.space;

import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.math.geometry.Cuboid;
import etomica.math.geometry.LineSegment;
import etomica.math.geometry.Polytope;
import etomica.math.geometry.Rectangle;
import etomica.math.geometry.Rectangular;
import etomica.simulation.Simulation;
import etomica.util.IRandom;

/**
 * Boundary that is in the shape of a rectangular parallelepiped.  
 * Periodicity in each direction is specified by subclass.
 */
public abstract class BoundaryRectangular extends Boundary implements BoundaryPeriodic {

    /**
     * Constructs cubic boundary of the given periodicity, using the space and default box-size
     * given by the Simulation. 
     */
    public BoundaryRectangular(Simulation sim, boolean[] periodicity) {
        this(sim.getSpace(), sim.getRandom(), periodicity, 30.0);
    }

    /**
     * Constructs cubic boundary of the given periodicity with each edge of length boxSize
     */
    public BoundaryRectangular(Space space, IRandom random, boolean[] periodicity, double boxSize) {
        this(space, random, periodicity, makeArray(space.D(), boxSize));
    }
    
    private static final double[] makeArray(int n, double d) {
        double[] array = new double[n];
        for (int i=0; i<n; i++) {
            array[i] = d;
        }
        return array;
    }

    /**
     * Constructs rectangular boundary of the given periodicity with edges given by the
     * values in the array boxSize.  Length of arrays must equal dimension of space.
     */
    public BoundaryRectangular(Space space, IRandom random, boolean[] periodicity, double[] boxSize) {
        super(space, makeShape(space));
        this.random = random;
        isPeriodic = (boolean[])periodicity.clone();
        dimensions = space.makeVector();
        dimensions.E(boxSize);
        
        temp = (IVectorRandom)space.makeVector();
        dimensionsCopy = space.makeVector();
        indexIterator = new IndexIteratorRectangular(space.D());
        needShift = new boolean[space.D()];//used by getOverflowShifts
        updateDimensions();
    }
    
    //used by constructors
    private static Polytope makeShape(Space space) {
        switch(space.D()) {
            case 1: return new LineSegment(space);
            case 2: return new Rectangle(space);
            case 3: return new Cuboid(space);
            default: throw new IllegalArgumentException("BoundaryRectangular not appropriate to given space");
        }
    }

    /**
     * Returns a vector with elements equal to the lengths of the edges of
     * the boundary.  The returned Vector does not represent the values internally,
     * so manipulation of the vector has no effect on this BoundaryRectangular instance.
     */
    public IVector getDimensions() {
        return dimensionsCopy;
    }
    
    public IVector getBoundingBox() {
        return dimensionsCopy;
    }

    /**
     * Returns a vector that describes a point selected uniformly within
     * the boundary.  The same Vector instance is returned with each call, with
     * a new random point each time.
     */
    public IVector randomPosition() {
        temp.setRandomCube(random);
        temp.TE(dimensions);
        return temp;
    }

    protected void updateDimensions() {
        dimensionsCopy.E(dimensions);
        ((Rectangular)shape).setEdgeLengths(dimensions);
    }

    /**
     * Sets the size and shape of the rectangular boundary.  Values are 
     * copied, so manipulation of the given vector has no subsequent effect
     * on this Boundary instance.
     */
    public void setDimensions(IVector v) {
        dimensions.E(v);
        updateDimensions();
    }

    /**
     * Returns the "volume" of the retangular region defined by this Boundary.
     * For a 2D and 1D spaces, this volume is actually an area and length, respectively.
     */
    public double volume() {
        return shape.getVolume();
    }

    /**
     * Returns the array that defines the directions of periodicity.  A true value
     * for each element indicates that the boundary is periodic in the correspoinding
     * direction.  The returned array is the internal representation, not a copy, so
     * changing it will affect the periodicity of the boundary.
     */
    public boolean[] getPeriodicity() {
        return isPeriodic;
    }

    public IVector[] getPeriodicVectors() {
        int n = 0;
        for (int i=0; i<isPeriodic.length; i++) {
            if (isPeriodic[i]) {
                n++;
            }
        }
        IVector[] vectors = new IVector[n];
        int d = 0;
        for  (int i=0; i<isPeriodic.length; i++) {
            if (isPeriodic[i]) {
                vectors[d] = space.makeVector();
                vectors[d].setX(i,dimensions.x(i));
                d++;
            }
        }
        return vectors;
    }
    public IndexIteratorSizable getIndexIterator() {
      int n = 0;
      for(int i=0; i<isPeriodic.length; i++)
        if(isPeriodic[i]) n++;
      return new IndexIteratorRectangular(n);
    }
    
    /**
     * Returns a set of image origins for a set of periodic image shells.  
     * The returned array is of dimension [(2*nShells+1)^D][D], where D
     * is the dimension of the space.
     */
    public double[][] imageOrigins(int nShells) {
        IVector workVector = space.makeVector();
        int shellFormula = (2 * nShells) + 1;
        int nImages = space.powerD(shellFormula) - 1;
        double[][] origins = new double[nImages][space.D()];
        indexIterator.setSize(shellFormula);
        indexIterator.reset();
        int k = 0;
        while(indexIterator.hasNext()) {
            int[] index = indexIterator.next();
            for (int i=0; i<space.D(); i++) {
                workVector.setX(i,index[i]);
            }
            workVector.PE(-(double)nShells);
            if(workVector.isZero()) continue;
            workVector.TE(dimensions);
            workVector.assignTo(origins[k++]);
        }
        return origins;
    }

    public float[][] getOverflowShifts(IVector rr, double distance) {
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
       double[] rrArray = new double[D];
       rr.assignTo(rrArray);
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
    
    private final IVectorRandom temp;
    protected final IVector dimensions;
    protected final IVector dimensionsCopy;
    private final IndexIteratorRectangular indexIterator;
    private final boolean[] needShift;
    protected boolean[] isPeriodic;
    protected final float[][] shift0 = new float[0][0];
    protected final IRandom random;

}
