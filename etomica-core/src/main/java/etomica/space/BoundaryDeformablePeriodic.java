/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.space;

import etomica.lattice.IndexIteratorRectangular;
import etomica.lattice.IndexIteratorSizable;
import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.Parallelogram;
import etomica.math.geometry.Parallelotope;
import etomica.math.geometry.Polytope;

/**
 * Boundary shaped as an arbitrary parallelepiped.  Applicable only for a 2D or 3D spaces.
 */

public class BoundaryDeformablePeriodic extends Boundary {

    /**
     * Make a cubic boundary with edges of length equal to the default boxSize and
     * periodic in every direction.
     */
	public BoundaryDeformablePeriodic(Space _space) {
		this(_space, 10.0);
	}

    /**
     * Make a cubic boundary of specified edge length and periodicity.
     */
	public BoundaryDeformablePeriodic(Space space, double boxSize) {
	    this(space, makeVectors(space, boxSize));
	}
	
    /**
     * Make a parallelepiped boundary of specified shape and periodicity.
     * 
     * @param space the governing space
     * @param vex array of vectors specifying directions for the parallelepiped edges
     * 
     *  @throws IllegalArgumentException if the dimension of space is not 2 or 3
     *  @throws IllegalArgumentException if the vex.length is not equal to the dimension of the space
     */
	public BoundaryDeformablePeriodic(Space space, Vector[] vex) {
        super(space, makeShape(space, vex));
        D = space.D();
        if(D != 2 && D != 3) {
            throw new IllegalArgumentException("BoundaryDeformablePeriodic is appropriate only for 2-D or 3-D spaces");
        }
        
        edgeVectors = new Vector[vex.length];
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i] = space.makeVector();
            edgeVectors[i].E(vex[i]);
        }

        h = space.makeTensor();
        hCopy = space.makeTensor();
        hInv = space.makeTensor();

        temp1 = space.makeVector();
        temp2 = space.makeVector();
        unit = space.makeVector();
        unit.E(1.0);
        half = space.makeVector();
        half.E(0.5);
        indexIterator = new IndexIteratorRectangular(space.D());
        periodicity = new boolean[space.D()];
        for (int i=0; i<periodicity.length; i++) {
            periodicity[i] = true;
        }
        update();
    }
    
    //used by constructor
    private static Vector[] makeVectors(Space space, double boxSize) {
        Vector[] vectors = new Vector[space.D()];
        for(int i=0; i<vectors.length; i++) {
            vectors[i] = space.makeVector();
            vectors[i].setX(i, boxSize);
        }
        return vectors;
    }
	
    //used only by constructor
    private static Polytope makeShape(Space space, Vector[] vex) {
        switch(space.D()) {
            case 2: return new Parallelogram(space, vex[0], vex[1]);
            case 3: return new Parallelepiped(space, vex[0], vex[1], vex[2]);
            default: throw new IllegalArgumentException("BoundaryDeformablePeriodic not appropriate to given space");
        }
    }

    /**
     * Returns the length of the sides of a rectangular box oriented in the lab
     * frame and in which the boundary is inscribed.
     */
	public Vector getBoxSize() {
        Vector[] vertices = shape.getVertices();
        temp1.E(vertices[0]);
        temp2.E(vertices[0]);
        for(int i=1; i<vertices.length; i++) {
            for (int j=0; j<space.D(); j++) {
                temp1.setX(j,Math.min(vertices[i].getX(j),temp1.getX(j)));
                temp2.setX(j,Math.max(vertices[i].getX(j),temp2.getX(j)));
            }
        }
        temp2.ME(temp1);
        return temp2;
    }

    /**
     * Returns the boundary tensor that fully describes the shape of the boundary.  Each
     * column of the tensor is an edge vector of the boundary.  Manipulation of the
     * returned tensor does not affect the state of the boundary. Shape-changing methods
     * (e.g., deform()) must be used to affect the boundary shape.
     */
    public Tensor getBoundaryTensor() {
        hCopy.E(h);
        return hCopy;
    }

    public Vector centralImage(Vector r) {
        temp1.E(r);
        for(int i=0; i<edgeVectors.length; i++) {
            temp1.PEa1Tv1(0.5,edgeVectors[i]);
        }
        hInv.transform(temp1);// position in terms of boundary-vector basis
        for (int i=0; i<space.D(); i++) {
            // remove any components that are likely nonzero due to roundoff
            if (Math.abs(temp1.getX(i)) < 1.e-10) {
                temp1.setX(i, 0);
            }
        }
        temp1.mod(unit);// shift to central image
        h.transform(temp1);//convert back to space-frame basis
        temp1.ME(r);//subtract r to return difference vector
        for(int i=0; i<edgeVectors.length; i++) {
            temp1.PEa1Tv1(-0.5,edgeVectors[i]);
        }
        return temp1;
    }

    public void nearestImage(Vector dr) {
        // To get out of the loop, we need to check n consecutive transformVectors
        // without applying any of them.  If we reach the end, then we wrap back around.
        for (int noTransformCount=0, i=0; noTransformCount<transformVectors.length; i++) {
            double dot = transformVectors[i].dot(dr)/tV2[i];
            if (Math.abs(dot) > halfTol) {
                dot = Math.round(dot);
                dr.PEa1Tv1(-dot, transformVectors[i]);
                // We transformed, but this now counts as a non-transfomrm -- we don't need to check it again.
                noTransformCount = 1;
            }
            else {
                noTransformCount++;
            }
            if (i==transformVectors.length-1) {
                // we finished a pass, but transformed a vector.
                // make another pass (we'll stop at the vector we transformed)
                i=-1;
            }
        }
        return;
    }


    /**
     * Applies the given deformation tensor to the boundary in its current
     * shape.
     * 
     * @throws IllegalArgumentException
     *             if the spatial dimension of the tensor is inconsistent with
     *             the dimension of the boundary
     */
    public void deform(etomica.space.Tensor deformationTensor) {
        if(deformationTensor.D() != D) {
            throw new IllegalArgumentException("Tensor dimension ("+deformationTensor.D()+") inconsistent with boundary dimension ("+D+")");           
        }
        for(int i=0; i<edgeVectors.length; i++) {
            deformationTensor.transform(edgeVectors[i]);
        }
        update();
    }
    
    /**
     * Sets the shape of the boundary as a rectangular parallelepiped (or
     * parallelogram) with edge lengths given by the elements of the given
     * vector.  To keep the edge lengths unchanged while deforming into
     * a rectangular shape, call this method with getDimensions() as its argument.
     * 
     * @throws IllegalArgumentException
     *             if the spatial dimension of the vector is inconsistent with
     *             the dimension of the boundary
     */
    public void setAsRectangular(Vector vector) {
        if(vector.getD() != D) {
            throw new IllegalArgumentException("Vector dimension ("+vector.getD()+") inconsistent with boundary dimension ("+D+")");
        }
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i].E(0.0);
            edgeVectors[i].setX(i, vector.getX(i));
        }
        update();
    }

    /**
     * Sets the shape of the boundary to a cube with edges equal to the given value.
     */
    public void setAsCube(double edgeLength) {
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i].E(0.0);
            edgeVectors[i].setX(i, edgeLength);
        }
        update();
    }
	/**
     * Scales each boundary edge so that its length equals the corresponding value
     * in the given vector.  Deformation of boundary is otherwise unchanged.
     */
	public void setBoxSize(Vector v) {
        if(!isPositive(v)) {
            throw new IllegalArgumentException("edge lengths must be greater than zero; attempt to set to "+v.toString());
        }
	    getBoxSize();
	    // temp2 is now current bounding box.
        temp1.E(v);
        temp1.DE(temp2);
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i].TE(temp1.getX(i));
        }
        update();
	}

	public void setEdgeVector(int i, Vector v){
	    edgeVectors[i].E(v);
        update();
	}
	
	
	/**
	 * Sets the boundary tensor to equal the given tensor.
	 */
	public void setDimensions(etomica.space.Tensor t) {
		t.assignTo(edgeVectors);
		update();
	}

	public double volume() {
		return volume;
	}

	public void setTruncationRadius(double newTruncationRadius) {
	    truncationRadius = newTruncationRadius;
	    update();
	}
    
	public double getTruncationRadius() {
	    return truncationRadius;
	}

    private boolean isPositive(Vector v) {
        for(int i=0; i<v.getD(); i++) {
            if(v.getX(i) <= 0.0) return false;
        }
        return true;
    }

    //method to update all auxiliary fields of tensor when edgeVectors are changed
    protected void update() {
//      h times a vector s gives coordinate r in lab frame; s elements are between 0,1 to be in box
        h.E(edgeVectors);//edge vectors in column format
        hInv.E(h);
        hInv.invert(); // to get point s in edgeVector frame, do hInv times r

        ((Parallelotope)shape).setEdgeVectors(edgeVectors);
        volume = shape.getVolume();

        transformVectors = new Vector[D];

        // add actual edges as transform vectors
        for (int i=0; i<D; i++) {
            transformVectors[i] = space.makeVector();
            transformVectors[i].E(edgeVectors[i]);
        }

        // test edge pairs
        for(int i=0; i<D-1; i++) {
            for (int k=i+1; k<D; k++) {
                double dot = edgeVectors[i].dot(edgeVectors[k]);
                if (Math.abs(dot) < 1e-10) continue;
                if(dot < -1e-10) {
                    temp1.Ev1Pv2(edgeVectors[i], edgeVectors[k]);
                }
                else if (dot > 1e-10) {
                    temp1.Ev1Mv2(edgeVectors[i], edgeVectors[k]);
                }
                testTransformVector(temp1);
            }
        }

        // test edge triplets 
        if (D>2) {
            double dot01 = edgeVectors[0].dot(edgeVectors[1]);
            double dot02 = edgeVectors[0].dot(edgeVectors[2]);
            double dot12 = edgeVectors[1].dot(edgeVectors[2]);

            double sum = 0.0;

            sum = dot01 + dot02 + dot12;//+0 +1 +2
            if(sum < -1e-10) {
                temp1.Ev1Pv2(edgeVectors[0], edgeVectors[1]);
                temp1.PE(edgeVectors[2]);
                testTransformVector(temp1);
            }

            sum = dot01 - dot02 - dot12;//+0 +1 -2
            if(sum < -1e-10) {
                temp1.Ev1Pv2(edgeVectors[0], edgeVectors[1]);
                temp1.ME(edgeVectors[2]);
                testTransformVector(temp1);
            }

            sum = -dot01 + dot02 - dot12;//+0 -1 +2
            if(sum < -1e-10) {
                temp1.Ev1Mv2(edgeVectors[0], edgeVectors[1]);
                temp1.PE(edgeVectors[2]);
                testTransformVector(temp1);
            }

            sum = -dot01 - dot02 + dot12;//+0 -1 -2
            if(sum < -1e-10) {
                temp1.Ev1Mv2(edgeVectors[0], edgeVectors[1]);
                temp1.ME(edgeVectors[2]);
                testTransformVector(temp1);
            }
        }

        if (tV2 == null || tV2.length != transformVectors.length) {
            tV2 = new double[transformVectors.length];
        }
        for (int i=0; i<transformVectors.length; i++) {
            tV2[i] = transformVectors[i].squared();
        }

        // we get called when the boundary changes, so fire inflate event now
        eventManager.inflate(this);
    }

    /**
     * We test the given vector (v) to see if 0.5*v would be transformed by any
     * of our existing transformVectors.  If so, then we don't need to keep v
     * as an additional transform vector.  We also check to see if vectors that
     * would be transformed could be within our truncation radius (if given).
     * If no vector transformed by only v would be within the truncation
     * radius, then v is not included as a new transform vector.
     */
    protected void testTransformVector(Vector v) {
        temp2.Ea1Tv1(0.5, v);
        for (int i=0; i<transformVectors.length; i++) {
            double dot = Math.abs(temp2.dot(transformVectors[i])/transformVectors[i].squared());
            if (dot > 1-halfTol) {
                return;
            }
        }

        if (truncationRadius < Double.POSITIVE_INFINITY) {
            // now check truncation distance
            temp2.Ea1Tv1(1-truncationRadius/Math.sqrt(v.squared()), v);
            // temp2 is now the periodic image of a vector on the edge of the
            // truncation radius.  If temp2 would be transformed by another
            // transformVector, then we don't need to keep v.
            for (int i=0; i<transformVectors.length; i++) {
                double dot = Math.abs(temp2.dot(transformVectors[i])/transformVectors[i].squared());
                if (dot > 1-halfTol) {
                    return;
                }
            }
        }

        // add v to our transform vectors
        transformVectors = (Vector[])etomica.util.Arrays.addObject(transformVectors, space.makeVector());
        transformVectors[transformVectors.length-1].E(v);
    }

    public Vector getEdgeVector(int d) {
        return edgeVectors[d];
    }
    
    public boolean getPeriodicity(int i) {
        return true;
    }
    
    public IndexIteratorSizable getIndexIterator() {
      return new IndexIteratorRectangular(edgeVectors.length);
    }
    
	public double[][] imageOrigins(int nShells) {
        int shellFormula = (2 * nShells) + 1;
        int nImages = space.powerD(shellFormula) - 1;
        if(nImages != origins.length) {
            origins = new double[nImages][space.D()];
        }
        indexIterator.setSize(shellFormula);
        indexIterator.reset();
        int k = 0;
        while (indexIterator.hasNext()) {
            int[] index = indexIterator.next();
            for (int i=0; i<index.length; i++) {
                temp1.setX(i,index[i]);
            }
            temp1.PE(-(double) nShells);
            if (temp1.isZero())
                continue;
            temp2.E(0.0);
            for (int i = 0; i < space.D(); i++) {
                temp2.PEa1Tv1(temp1.getX(i), edgeVectors[i]);
            }
            temp2.assignTo(origins[k++]);
        }
        return origins;
    }//end of imageOrigins

    int shiftX, shiftY, shiftZ;
    private double volume;
    private Tensor h;
    private Tensor hCopy;
    private final Tensor hInv;
    private final Vector[] edgeVectors;
    private final Vector temp1;
    private final Vector temp2;
    private final boolean[] periodicity;
    private final Vector unit;
    private final Vector half;
    private final int D;
    private final IndexIteratorRectangular indexIterator;
    private double[][] origins = new double[0][];
    private final static double halfTol = 0.50000000001;
    protected Vector[] transformVectors;
    protected double[] tV2;
    protected double truncationRadius = Double.POSITIVE_INFINITY;

    private static final long serialVersionUID = 1L;
}
