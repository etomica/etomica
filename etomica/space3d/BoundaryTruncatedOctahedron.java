package etomica.space3d;

import etomica.exception.MethodNotImplementedException;
import etomica.lattice.IndexIteratorSequential;
import etomica.math.geometry.Plane;
import etomica.math.geometry.Polygon;
import etomica.math.geometry.Polyhedron;
import etomica.math.geometry.TruncatedOctahedron;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryPeriodic;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * This class enables creation of a periodic truncated-octahedron boundary.
 * Truncated octahedrons are spacefillers; unlike spheres, about which voids
 * form when placed in an array, truncated octahedrons can be flush with other
 * truncated octahedrons at every face when grouped in an array. Moreover,
 * truncated octahedrons are considerably more sphere-like than cubes, and using
 * a truncated octahedron as a boundary greatly increases the isotropy of the
 * system.
 * <p>
 * There are no subclasses for periodic and nonperiodic boundary conditions for
 * this class because nonperiodic boundary conditions are not applicable for
 * truncated octahedrons.
 */
/*
 * History Created on Jan 24, 2005 by kofke
 */
public class BoundaryTruncatedOctahedron extends Boundary implements
        BoundaryPeriodic {

    public BoundaryTruncatedOctahedron(Simulation sim) {
        this(sim.space, sim.getDefaults().boxSize);
    }
    public BoundaryTruncatedOctahedron(Space space, double boxSize) {
        super(space, new TruncatedOctahedron(space));
        isPeriodic = new boolean[space.D()];
        for (int i = 0; i < space.D(); i++)
            isPeriodic[i] = true;
        dimensions = space.makeVector();
        dimensions.E(boxSize);
        rrounded = space.makeVector();
        intoTruncatedOctahedron = space.makeVector();
        intoContainingCube = space.makeVector();
        dimensionsCopy = space.makeVector();
        dimensionsHalf = space.makeVector();
        indexIterator = new IndexIteratorSequential(space.D());
        needShift = new boolean[space.D()];//used by getOverflowShifts
        updateDimensions();
    }

    public boolean[] getPeriodicity() {
        return isPeriodic;
    }

    public final etomica.space.Vector dimensions() {
        return dimensionsCopy;
    }

    public etomica.space.Vector randomPosition() {
        throw new MethodNotImplementedException();
        //temp.setRandom(dimensions);
        //return temp;
    }

    private final void updateDimensions() {
        dimensionsHalf.Ea1Tv1(0.5, dimensions);
        dimensionsCopy.E(dimensions);
        ((TruncatedOctahedron) shape).setContainingCubeEdgeLength(dimensions
                .x(0));
    }

    public void setDimensions(etomica.space.Vector v) {
        dimensions.E(v);
        updateDimensions();
    }

    public double[][] imageOrigins(int nShells) {
        if(nShells == 0) {
            return new double[0][];
        } else if(nShells == 1) {
            Polygon[] faces = ((Polyhedron)shape).getFaces();
            double[][] origins = new double[faces.length][];
            double multiplier = ((TruncatedOctahedron)shape).getContainingCubeEdgeLength();
            for(int i=0; i<faces.length; i++) {
                Vector[] vertices = faces[i].getVertices();
                plane.setThreePoints((Vector3D)vertices[0], (Vector3D)vertices[1], (Vector3D)vertices[2]);
                plane.getNormalVector(normal);
                normal.TE(multiplier);
                origins[i] = normal.toArray();
            }
            return origins;
        }
        //algorithm for nShells > 1 misses many of the images (those through the hexagon faces)
        Vector workVector = space.makeVector();
        int shellFormula = (2 * nShells) + 1;
        int nImages = space.powerD(shellFormula) - 1;
        double[][] origins = new double[nImages][space.D()];
        indexIterator.setSize(shellFormula);
        indexIterator.reset();
        int k = 0;
        while (indexIterator.hasNext()) {
            int[] index = indexIterator.next();
            workVector.E(index);
            workVector.PE(-(double) nShells);
            if (workVector.isZero())
                continue;
            workVector.TE(dimensions);
            workVector.assignTo(origins[k++]);
        }
        return origins;
    }

    public float[][] getOverflowShifts(Vector rr, double distance) {
        int D = space.D();
        int numVectors = 1;
        for (int i = 1; i < D; i++) {
            if ((rr.x(i) - distance < 0.0)
                    || (rr.x(i) + distance > dimensions.x(i))) {
                //each previous vector will need an additional copy in this
                // dimension
                numVectors *= 2;
                //remember that
                needShift[i] = true;
            } else {
                needShift[i] = false;
            }
        }

        if (numVectors == 1)
            return shift0;

        float[][] shifts = new float[numVectors][D];
        double[] rrArray = rr.toArray();
        for (int i = 0; i < D; i++) {
            shifts[0][i] = (float) rrArray[i];
        }
        int iVector = 1;

        for (int i = 0; iVector < numVectors; i++) {
            if (!needShift[i]) {
                //no shift needed for this dimension
                continue;
            }
            double delta = -dimensions.x(i);
            if (rr.x(i) - distance < 0.0) {
                delta = -delta;
            }
            //copy all previous vectors and apply a shift of delta to the
            // copies
            for (int jVector = 0; jVector < iVector; jVector++) {
                for (int j = 0; j < D; j++) {
                    shifts[jVector + iVector][j] = shifts[jVector][j];
                }
                shifts[jVector + iVector][i] += delta;
            }
            iVector *= 2;
        }
        return shifts;
    }

    public void nearestImage(Vector dr) {
        dr.PEa1Tv1(0.5, dimensions);
        dr.PE(centralImage(dr));
        dr.PEa1Tv1(-0.5, dimensions);

        //        double n =
        // ((TruncatedOctahedron)shape).getContainingCubeEdgeLength();
        //        intoContainingCube.EModShift(dr, dimensions);
        //        rrounded.Ev1Pv2(dr,intoContainingCube);
        //        rrounded.TE(1./n);
        //        int aint =
        // (int)(4.0/3.0*(Math.abs(rrounded.x(0))+Math.abs(rrounded.x(1))+Math.abs(rrounded.x(2))));
        //        double corr = 0.5 * n * aint;
        //        if(corr != 0.0) {
        //            if(rrounded.x(0) > 0) intoTruncatedOctahedron.setX(0,
        // intoContainingCube.x(0) - corr);
        //                else intoTruncatedOctahedron.setX(0, intoContainingCube.x(0) + corr);
        //            if(rrounded.x(1) > 0) intoTruncatedOctahedron.setX(1,
        // intoContainingCube.x(1) - corr);
        //                else intoTruncatedOctahedron.setX(1, intoContainingCube.x(1) + corr);
        //            if(rrounded.x(2) > 0) intoTruncatedOctahedron.setX(2,
        // intoContainingCube.x(2) - corr);
        //                else intoTruncatedOctahedron.setX(2, intoContainingCube.x(2) + corr);
        //        } else {
        //            intoTruncatedOctahedron.E(intoContainingCube);
        //        }
        //        dr.PE(intoTruncatedOctahedron);
    }

    public Vector centralImage(Vector r) {
        double n = ((TruncatedOctahedron) shape).getContainingCubeEdgeLength();
        intoContainingCube.EModShift(r, dimensions);
        rrounded.Ev1Pv2(r, intoContainingCube);
        rrounded.TE(1. / n);
        rrounded.PE(-0.5);
        int aint = (int) (4.0 / 3.0 * (Math.abs(rrounded.x(0))
                + Math.abs(rrounded.x(1)) + Math.abs(rrounded.x(2))));
        double corr = 0.5 * n * aint;
        intoTruncatedOctahedron.E(intoContainingCube);

        if (corr != 0.0) {
            intoTruncatedOctahedron.PEa1SGNv1(-corr, rrounded);
        }
        return intoTruncatedOctahedron;
    }

    protected final Vector intoTruncatedOctahedron;
    protected final Vector rrounded;
    protected final Vector intoContainingCube;
    protected final Vector dimensions;
    protected final Vector dimensionsCopy;
    protected final Vector dimensionsHalf;
    private final IndexIteratorSequential indexIterator;
    private final boolean[] needShift;
    protected final boolean[] isPeriodic;
    protected final float[][] shift0 = new float[0][0];
    protected float[][] shift;
    private final Plane plane = new Plane();
    private final Vector3D normal = new Vector3D();

}
