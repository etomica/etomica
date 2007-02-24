package etomica.space;

import java.io.Serializable;

import etomica.lattice.IndexIteratorSequential;
import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.Parallelogram;
import etomica.math.geometry.Parallelotope;
import etomica.math.geometry.Polytope;
import etomica.simulation.Simulation;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;

/**
 * Boundary shaped as an arbitrary parallelepiped.  Applicable only for a 2D or 3D spaces.
 */

//nan needs a cleanup of getOverflowShifts.
public class BoundaryDeformablePeriodic extends Boundary {

    /**
     * Make a cubic boundary with edges of length equal to the default boxSize and
     * periodic in every direction.
     */
	public BoundaryDeformablePeriodic(Simulation sim) {
		this(sim.getSpace(), sim.getDefaults().boxSize);
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
	public BoundaryDeformablePeriodic(Space space,  Vector[] vex) {
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

        if(D == 2) { 
            edgePairTransforms = new PeriodicTransform2[1];
            edgeTripletTransforms = new PeriodicTransform3[0];
            edgePairTransforms[0] = new PeriodicTransform2(edgeVectors, 0, 1);
        } else {//D == 3
            edgePairTransforms = new PeriodicTransform2[3];
            edgeTripletTransforms = new PeriodicTransform3[2];
            edgePairTransforms[0] = new PeriodicTransform2(edgeVectors, 1, 2);
            edgePairTransforms[1] = new PeriodicTransform2(edgeVectors, 0, 2);
            edgePairTransforms[2] = new PeriodicTransform2(edgeVectors, 0, 1);
            edgeTripletTransforms[0] = new PeriodicTransform3(edgeVectors, true);
            edgeTripletTransforms[1] = new PeriodicTransform3(edgeVectors, false);
            
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
        dimensions = space.makeVector();
        dimensionsCopy = space.makeVector();
        dimensionsHalf = space.makeVector();
        indexIterator = new IndexIteratorSequential(space.D());
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
            case 2: return new Parallelogram(space, (Vector2D)vex[0], (Vector2D)vex[1]);
            case 3: return new Parallelepiped(space, (Vector3D)vex[0], (Vector3D)vex[1], (Vector3D)vex[2]);
            default: throw new IllegalArgumentException("BoundaryDeformablePeriodic not appropriate to given space");
        }
    }

    /**
     * Returns a vector with each element equal to the length of the corresponding
     * edge of the boundary.
     */
	public Vector getDimensions() {
        dimensionsCopy.E(dimensions);
        return dimensionsCopy;
    }

    /**
     * Returns the boundary tensor that fully describes the shape of the boundary.  Each
     * column of the tensor is an edge vector of the boundary.  Manipulation of the
     * returned tensor does not affect the state of the boundary. Shape-changing methods
     * (e.g., deform()) must be used to affect the boundary shape.
     */
    public etomica.space.Tensor boundaryTensor() {
        hCopy.E(h);
        return hCopy;
    }

    /**
     * Returns a point selected randomly from within the boundary region.
     */
    public etomica.space.Vector randomPosition() {
        temp1.setRandom(1.0);
        temp1.transform(h);
        return temp1;
    }

    public Vector getBoundingBox() {
        Vector[] vertices = shape.getVertices();
        temp1.E(vertices[0]);
        temp2.E(vertices[0]);
        for(int i=1; i<vertices.length; i++) {
            temp1.minE(vertices[i]);
            temp2.maxE(vertices[i]);
        }
        temp2.ME(temp1);
        return temp2;
    }
    
    public Vector centralImage(Vector r) {
        temp1.E(r);
        for(int i=0; i<edgeVectors.length; i++) {
            temp1.PEa1Tv1(0.5,edgeVectors[i]);
        }
        temp1.transform(hInv);// position in terms of boundary-vector basis
        temp1.truncate(1.e-10);// remove any components that are likely nonzero due to roundoff
        temp1.mod(unit);// shift to central image
        temp1.transform(h);//convert back to space-frame basis
        temp1.ME(r);//subtract r to return difference vector
        for(int i=0; i<edgeVectors.length; i++) {
            temp1.PEa1Tv1(-0.5,edgeVectors[i]);
        }
        return temp1;
    }

    //needs work to improve efficiency; may be incorrect for extremely deformed boundaries
	public void nearestImage(Vector dr) { 
	    
        boolean transformed = false;
       // temp1.E(dr);
       // temp1.transform(hInv);//transform to edge-vector basis
//        for(int i=0; i<edgePairTransforms.length; i++) {
//            transformed = edgePairTransforms[i].transform(temp1);
//            if(transformed) {
//                temp1.transform(h);//convert back to space-frame basis
//                dr.E(temp1);
//                return;
//            }
//        }
//
        //transform into reciprocal lattice unit cell
        double dot = 0.0;
        //temp1.E(dr);
        //temp1.transform(hInv);
        //System.out.println(temp1.toString());
        do {
            transformed = false;
            for(int i=0; i<edgeVectors.length; i++) {
                dot = dr.dot(edgeVectors[i])/edgeVectors[i].squared();
                //System.out.println("dot "+dot);
                if(dot > 0.5) {
                    //System.out.println("subtract "+edgeVectors[i]);
                    do {
                        dr.ME(edgeVectors[i]);
                        dot -= 1.0;
                        //System.out.println("minus "+edgeVectors[i]);
                    } while(dot > 0.5);
                    transformed = true;
                } else if(dot < -0.5) {
                    //System.out.println("add "+edgeVectors[i]);
                    do {
                        dr.PE(edgeVectors[i]);
                        dot += 1.0;
                        //System.out.println("plus"+edgeVectors[i]);
                    } while(dot < -0.5);
                    transformed = true;
                }
                //dot = dr.dot(edgeVectors[i])/edgeVectors[i].squared();
    //            if(dot > 0.5 || dot < -0.5) {
    //                System.out.println("busted");
    //            }
            }
        } while(transformed);
        //temp1.E(dr);
        //temp1.transform(hPrimeInv);// position in terms of reciprocal-vector basis
        //temp1.PE(half);
        //temp1.mod(unit);
        //temp1.ME(half);
        
        transformed = false;
        for(int i=0; i<edgeTripletTransforms.length; i++) {
            transformed = edgeTripletTransforms[i].transform(dr);
        }
//        if(transformed) {
//            return;
//        }
        
        //temp1.transform(hPrime);//transform back to lab frame
        //temp1.transform(hInv);//transform to edge-vector basis
        
        //System.out.println(temp1.toString());
        
        do {
            transformed = false;
            for(int i=0; i<edgePairTransforms.length; i++) {
                boolean trans = edgePairTransforms[i].transform(dr);
                transformed |= trans;
                if(trans) {
                    dot = dr.dot(edgeVectors[i])/edgeVectors[i].squared();
                    //System.out.println("dot "+dot);
                    if(dot > 0.5) {
                        //System.out.println("subtract "+edgeVectors[i]);
                            dr.ME(edgeVectors[i]);
                            //System.out.println("minus "+edgeVectors[i]);
                    } else if(dot < -0.5) {
                        //System.out.println("add "+edgeVectors[i]);
                            dr.PE(edgeVectors[i]);
                            //System.out.println("plus"+edgeVectors[i]);
                    }
                }
            }
        } while(transformed);
        do {
            transformed = false;
            for(int i=0; i<edgeVectors.length; i++) {
                dot = dr.dot(edgeVectors[i])/edgeVectors[i].squared();
                //System.out.println("dot "+dot);
                if(dot > 0.5) {
                    //System.out.println("subtract "+edgeVectors[i]);
                    do {
                        dr.ME(edgeVectors[i]);
                        dot -= 1.0;
                        //System.out.println("minus "+edgeVectors[i]);
                    } while(dot > 0.5);
                    transformed = true;
                } else if(dot < -0.5) {
                    //System.out.println("add "+edgeVectors[i]);
                    do {
                        dr.PE(edgeVectors[i]);
                        dot += 1.0;
                        //System.out.println("plus"+edgeVectors[i]);
                    } while(dot < -0.5);
                    transformed = true;
                }
                //dot = dr.dot(edgeVectors[i])/edgeVectors[i].squared();
    //            if(dot > 0.5 || dot < -0.5) {
    //                System.out.println("busted");
    //            }
            }
        } while(transformed);
        
        transformed = false;
        for(int i=0; i<edgeTripletTransforms.length; i++) {
            transformed = edgeTripletTransforms[i].transform(dr);
        }
        if(transformed) {
            return;
        }
        
        //System.out.println();

       //System.out.println(temp1.toString());
        //temp1.transform(h);//convert back to space-frame basis
        //dr.E(temp1);

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
            edgeVectors[i].transform(deformationTensor);
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
        if(vector.D() != D) {
            throw new IllegalArgumentException("Vector dimension ("+vector.D()+") inconsistent with boundary dimension ("+D+")");
        }
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i].E(0.0);
            edgeVectors[i].setX(i, vector.x(i));
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
	public void setDimensions(etomica.space.Vector v) {
        if(!isPositive(v)) {
            throw new IllegalArgumentException("edge lengths must be greater than zero; attempt to set to "+v.toString());
        }
        temp1.E(v);
        temp1.DE(dimensions);
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i].TE(temp1.x(i));
        }
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
    
    private boolean isPositive(Vector v) {
        for(int i=0; i<v.D(); i++) {
            if(v.x(i) <= 0.0) return false;
        }
        return true;
    }

    //method to update all auxiliary fields of tensor when edgeVectors are changed
    private void update() {
//      h times a vector s gives coordinate r in lab frame; s elements are between 0,1 to be in box
        h.E(edgeVectors);//edge vectors in column format
        hInv.E(h);
        hInv.inverse(); // to get point s in edgeVector frame, do hInv times r
        
        for(int i=0; i<edgeVectors.length; i++) {
            dimensions.setX(i, Math.sqrt(edgeVectors[i].squared()));
        }
        dimensionsHalf.Ea1Tv1(0.5, dimensions);
        ((Parallelotope)shape).setEdgeVectors(edgeVectors);
        volume = shape.getVolume();
        
        for(int i=0; i<edgePairTransforms.length; i++) {
            edgePairTransforms[i].update();
        }
        for(int i=0; i<edgeTripletTransforms.length; i++) {
            edgeTripletTransforms[i].update();
        }
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
            temp1.E(index);
            temp1.PE(-(double) nShells);
            if (temp1.isZero())
                continue;
            temp2.E(0.0);
            for (int i = 0; i < space.D(); i++) {
                temp2.PEa1Tv1(temp1.x(i), edgeVectors[i]);
            }
            temp2.assignTo(origins[k++]);
        }
        return origins;
    }//end of imageOrigins

	public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {
		throw new RuntimeException(
				"BoundaryDeformablePeriodic.getOverflowShifts not implmented");
		/*
		 * shiftX = 0; shiftY = 0; shiftZ = 0; r = (Vector)rr;
		 * 
		 * if(r.x-distance < 0.0) {shiftX = +1;} else if(r.x+distance >
		 * dimensions.x) {shiftX = -1;}
		 * 
		 * if(r.y-distance < 0.0) {shiftY = +1;} else if(r.y+distance >
		 * dimensions.y) {shiftY = -1;}
		 * 
		 * if(r.z-distance < 0.0) {shiftZ = +1;} else if(r.z+distance >
		 * dimensions.z) {shiftZ = -1;}
		 * 
		 * if((shiftX == 0) && (shiftY == 0) && (shiftZ == 0)) { shift = shift0; }
		 * else if((shiftX != 0) && (shiftY == 0) && (shiftZ == 0)) { shift =
		 * new float[1][D]; shift[0][0] = (float)(shiftX*dimensions.x); } else
		 * if((shiftX == 0) && (shiftY != 0) && (shiftZ == 0)) { shift = new
		 * float[1][D]; shift[0][1] = (float)(shiftY*dimensions.y); } else
		 * if((shiftX == 0) && (shiftY == 0) && (shiftZ != 0)) { shift = new
		 * float[1][D]; shift[0][2] = (float)(shiftZ*dimensions.z); } else
		 * if((shiftX != 0) && (shiftY != 0) && (shiftZ == 0)) { shift = new
		 * float[3][D]; shift[0][0] = (float)(shiftX*dimensions.x); shift[1][1] =
		 * (float)(shiftY*dimensions.y); shift[2][0] = shift[0][0]; shift[2][1] =
		 * shift[1][1]; } else if((shiftX != 0) && (shiftY == 0) && (shiftZ !=
		 * 0)) { shift = new float[3][D]; shift[0][0] =
		 * (float)(shiftX*dimensions.x); shift[1][2] =
		 * (float)(shiftZ*dimensions.z); shift[2][0] = shift[0][0]; shift[2][2] =
		 * shift[1][2]; } else if((shiftX == 0) && (shiftY != 0) && (shiftZ !=
		 * 0)) { shift = new float[3][D]; shift[0][1] =
		 * (float)(shiftY*dimensions.y); shift[1][2] =
		 * (float)(shiftZ*dimensions.z); shift[2][1] = shift[0][1]; shift[2][2] =
		 * shift[1][2]; } else if((shiftX != 0) && (shiftY != 0) && (shiftZ !=
		 * 0)) { shift = new float[7][D]; shift[0][0] =
		 * (float)(shiftX*dimensions.x); shift[1][1] =
		 * (float)(shiftY*dimensions.y); shift[2][2] =
		 * (float)(shiftZ*dimensions.z); shift[3][0] = shift[0][0]; shift[3][1] =
		 * shift[1][1]; shift[4][1] = shift[1][1]; shift[4][2] = shift[2][2];
		 * shift[5][0] = shift[0][0]; shift[5][2] = shift[2][2]; shift[6][0] =
		 * shift[0][0]; shift[6][1] = shift[1][1]; shift[6][2] = shift[2][2]; }
		 * 
		 * return(shift);
		 */
	}//end of getOverflowShifts

    int shiftX, shiftY, shiftZ;
    private double volume;
    private Tensor h;
    private Tensor hCopy;
    private final Tensor hInv;
    private final Vector[] edgeVectors;
    private final Vector temp1, temp2;
    private final Vector dimensions;
    private final Vector dimensionsCopy;
    private final Vector dimensionsHalf;
    private final Vector unit;
    private final Vector half;
    private final int D;
    private final PeriodicTransform2[] edgePairTransforms;
    private final PeriodicTransform3[] edgeTripletTransforms;
    private final IndexIteratorSequential indexIterator;
    private double[][] origins = new double[0][];

    private static final long serialVersionUID = 1L;
    
    private static abstract class PeriodicTransform implements Serializable {
        
        PeriodicTransform(int D) {
            transformVector = Space.makeVector(D);
        }
        
        boolean transform(Vector dr) {
            double dot = dr.dot(transformVector)/t2;
            if(dot > 0.5) {
                do {
                    //System.out.println(n+"minus"+dr+transformVector.toString());
                    dr.ME(transformVector);
                    dot -= 1;
                } while(dot > 0.5); 
                return true;
            } else if(dot < -0.5) {
                do {
                    //System.out.println(n+"plus"+dr+transformVector.toString());
                    dr.PE(transformVector);
                    dot += 1;
                } while(dot < -0.5);
                return true;
            }
            return false;
        }
        
        
        abstract void update();
        
        protected final Vector transformVector;
        protected double t2;
        protected int n;
   }
    
    private static class PeriodicTransform2 extends PeriodicTransform {
        
        PeriodicTransform2(Vector[] edgeVectors, int k0, int k1) {
            super(edgeVectors.length);
            edge0 = edgeVectors[k0];
            edge1 = edgeVectors[k1];
            n = 2;
        }
        
        void update() {
            if(edge0.dot(edge1) < 0) {
                transformVector.Ev1Pv2(edge0,edge1);
            } else {
                transformVector.Ev1Mv2(edge0,edge1);
            }
            t2 = transformVector.squared();
        }
        
        private final Vector edge0, edge1;
        private static final long serialVersionUID = 1L;
    }

    private static class PeriodicTransform3 extends PeriodicTransform {
        
        PeriodicTransform3(Vector[] edgeVectors, boolean getFirst) {
            super(edgeVectors.length);
            edge0 = edgeVectors[0];
            edge1 = edgeVectors[1];
            edge2 = edgeVectors[2];
            this.getFirst = getFirst;
            n = 3;
        }
        
        void update() {
            boolean getNext = getFirst;
            
            double dot01 = edge0.dot(edge1);
            double dot02 = edge0.dot(edge2);
            double dot12 = edge1.dot(edge2);
            
            double sum = 0.0;
            
            sum = dot01 + dot02 + dot12;//+0 +1 +2
            if(sum < 0) {
                if(getNext) {
                    transformVector.Ev1Pv2(edge0, edge1);
                    transformVector.PE(edge2);
                    t2 = transformVector.squared();
                    return;
                }
                getNext = true;
            }
            
            sum = dot01 - dot02 - dot12;//+0 +1 -2
            if(sum < 0) {
                if(getNext) {
                    transformVector.Ev1Pv2(edge0, edge1);
                    transformVector.ME(edge2);
                    t2 = transformVector.squared();
                    return;
                }
                getNext = true;
            }
            
            sum = -dot01 + dot02 - dot12;//+0 -1 +2
            if(sum < 0) {
                if(getNext) {
                    transformVector.Ev1Mv2(edge0, edge1);
                    transformVector.PE(edge2);
                    t2 = transformVector.squared();
                    return;
                }
                getNext = true;
            }
            
            sum = -dot01 - dot02 + dot12;//+0 -1 -2
            if(sum < 0) {
                if(getNext) {
                    transformVector.Ev1Mv2(edge0, edge1);
                    transformVector.ME(edge2);
                    t2 = transformVector.squared();
                    return;
                }
            }
            if(sum != 0) throw new RuntimeException("i really should't be here");
        }
        
        private final Vector edge0, edge1, edge2;
        private final boolean getFirst;
        private static final long serialVersionUID = 1L;
    }

}//end of BoundaryDeformablePeriodic
