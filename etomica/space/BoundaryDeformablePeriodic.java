package etomica.space;

import etomica.lattice.IndexIteratorSequential;
import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.Parallelogram;
import etomica.math.geometry.Parallelotope;
import etomica.math.geometry.Polytope;
import etomica.simulation.Simulation;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;

/**
 * Boundary shaped as an arbitrary parallelepiped.  Applicable only for a 3D space.
 */

//nan needs a cleanup of imageOrigins & getOverflowShifts.
/**
 * Warning, this class assumes a rectangular system!!
 * @author nancycribbin
 *
 */
public class BoundaryDeformablePeriodic extends Boundary implements BoundaryPeriodic {


    /**
     * Make a cubic boundary with edges of length equal to the default boxSize and
     * periodic in every direction.
     */
	public BoundaryDeformablePeriodic(Simulation sim) {
		this(sim.space, makePeriodicity(sim.space.D()), sim.getDefaults().boxSize);
	}

    /**
     * Make a cubic boundary of specified edge length and periodicity.
     */
	public BoundaryDeformablePeriodic(Space space, boolean[] periodicity, double boxSize) {
	    this(space, periodicity, makeVectors(space, boxSize));	
	}
	
    /**
     * Make a parallelepiped boundary of specified shape and periodicity.
     * 
     * @param space the governing space (must be Space3D)
     * @param periodicity array indicating for each direction whether it is periodic (true) or not (false)
     * @param vex array of three vectors specifying directions for three of the parallelepiped edges 
     */
	public BoundaryDeformablePeriodic(Space space, boolean[] periodicity, Vector[] vex) {
        super(space, makeShape(space, vex));
        edgeVectors = (Vector[])vex.clone();
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i] = (Vector)vex[i].clone();
        }
        isPeriodic = (boolean[]) periodicity.clone();

        h = space.makeTensor();
        hCopy = space.makeTensor();
        hInv = space.makeTensor();

        temp1 = space.makeVector();
        temp2 = space.makeVector();
        unit = space.makeVector();
        unit.E(1.0);
        dimensions = space.makeVector();
        dimensionsCopy = space.makeVector();
        dimensionsHalf = space.makeVector();
        tempVector = space.makeVector();
        nearestDr = space.makeVector();
        tempTensor = space.makeTensor();
        indexIterator = new IndexIteratorSequential(space.D());
//        tempVector.E(vex[0].P(vex[1].P(vex[2])));
//        setDimensions(tempVector);
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

    //used only by constructor
	private static boolean[] makePeriodicity(int D) {
		boolean[] isPeriodic = new boolean[D];
		for (int i = 0; i < D; i++) {
			isPeriodic[i] = true;
		}
		return isPeriodic;
	}

    
    public boolean[] getPeriodicity() {
        return isPeriodic;
    }

	public Vector getDimensions() {
        dimensionsCopy.E(dimensions);
        return dimensionsCopy;
    }

    public etomica.space.Tensor boundaryTensor() {
        hCopy.E(h);
        return hCopy;
    }

    public etomica.space.Vector randomPosition() {
        temp1.setRandom(1.0);
        temp1.transform(h);
        return temp1;
    }

    private void update() {
        h.E(edgeVectors);
        hInv.E(h);
        hInv.inverse(); // to get s, inverse of h times r, which is a position
        for(int i=0; i<edgeVectors.length; i++) {
            dimensions.setX(i, Math.sqrt(edgeVectors[i].squared()));
        }
        dimensionsHalf.Ea1Tv1(0.5, dimensions);
        ((Parallelotope)shape).setEdgeVectors(edgeVectors);
    }

	public Vector centralImage(Vector r) {
        temp1.E(r);
        temp1.transform(hInv);// position in terms of boundary-vector basis
        temp1.truncate(1.e-10);// remove any components that are likely nonzero due to roundoff
        temp1.mod(unit);// shift to central image
        temp1.transform(h);//convert back to space-frame basis
        temp1.ME(r);//subtract r to return difference vector
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
    
    public Vector getCenter() {
        temp1.E(edgeVectors[0]);
        for(int i=1; i<edgeVectors.length; i++) {
            temp1.PE(edgeVectors[i]);
        }
        temp1.TE(0.5);
        return temp1;
    }

	public void setUseBruteForceNearestImageMethod(boolean bb) {
		useBruteForceNearestImageMethod = bb;
	}
//skkwak next!!!
	
	/**
	 * Checks for the nearest image.  The argument vector is changed by this method.
	 */
	public void nearestImage(Vector dr2) { // dr2 is not dr^2
	    //the stuff after the first if is the "new technique"
		if (!useBruteForceNearestImageMethod) {
			if (!constantStrain) {
				update();
			} //in Main Class, makeStrainTensor() must be called then set up
			  // boolean of contantStrain!!!

			double xi = dr2.x(0) * hInv.component(0, 0) + dr2.x(1)
					* hInv.component(0, 1) + dr2.x(2) * hInv.component(0, 2);
			double eta = dr2.x(0) * hInv.component(1, 0) + dr2.x(1)
					* hInv.component(1, 1) + dr2.x(2) * hInv.component(1, 2);
			double zeta = dr2.x(0) * hInv.component(2, 0) + dr2.x(1)
					* hInv.component(2, 1) + dr2.x(2) * hInv.component(2, 2);

			//This series of if statements sets xi, eta, and zeta equal to 0 if they are close enough to it (i.e. very very small).
			if (xi > -1.0e-10 && xi < 1.0e-10) {
				xi = 0;
			}
			if (eta > -1.0e-10 && eta < 1.0e-10) {
				eta = 0;
			}
			if (zeta > -1.0e-10 && zeta < 1.0e-10) {
				zeta = 0;
			}

			//Increment xi, eta, & zeta
			xi += 0.5;
			eta += 0.5;
			zeta += 0.5;

			//Move xi, eta, and zeta between 0 and 1, if they are not already 
			while (xi < 0) {
				xi += 1;
			}
			while (xi > 1) {
				xi -= 1;
			}
			while (eta < 0) {
				eta += 1;
			}
			while (eta > 1) {
				eta -= 1;
			}
			while (zeta < 0) {
				zeta += 1;
			}
			while (zeta > 1) {
				zeta -= 1;
			}

			//Decrement xi, eta, & zeta
			xi -= 0.5;
			eta -= 0.5;
			zeta -= 0.5;

			dr2.setX(0, xi * h.component(0, 0) + eta
					* h.component(1, 0) + zeta
					* h.component(2, 0));
			dr2.setX(1, xi * h.component(0, 1) + eta
					* h.component(1, 1) + zeta
					* h.component(2, 1));
			dr2.setX(2, xi * h.component(0, 2) + eta
					* h.component(1, 2) + zeta
					* h.component(2, 2));

			//                    } //inner if
		} else {
			if (!constantStrain) {
				update();
			} //in Main Class, makeStrainTensor() must be called then set up
			  // boolean of contantStrain!!!

			oldRefDistance = Math.sqrt(dr2.squared());
			nearestDr.E(dr2);

			double xi = dr2.x(0) * hInv.component(0, 0) + dr2.x(1)
					* hInv.component(0, 1) + dr2.x(2) * hInv.component(0, 2);
			double eta = dr2.x(0) * hInv.component(1, 0) + dr2.x(1)
					* hInv.component(1, 1) + dr2.x(2) * hInv.component(1, 2);
			double zeta = dr2.x(0) * hInv.component(2, 0) + dr2.x(1)
					* hInv.component(2, 1) + dr2.x(2) * hInv.component(2, 2);

			workXEZ[0] = xi;
			workXEZ[1] = eta;
			workXEZ[2] = zeta;

			if (xi > 0) {
				workXEZ[0] -= 1;
			} else {
				workXEZ[0] += 1;
			}
			if (eta > 0) {
				workXEZ[1] -= 1;
			} else {
				workXEZ[1] += 1;
			}
			if (zeta > 0) {
				workXEZ[2] -= 1;
			} else {
				workXEZ[2] += 1;
			}

			dr2.E(compareDistance(xi, eta, zeta, workXEZ));

		}//end of if and else
	}

	public etomica.space.Vector compareDistance(double xi, double eta,
			double zeta, double[] xez) {
		double X = xi, E = eta, Z = zeta;
		for (int i = 0; i < 7; i++) {
			if (i == 0) {
				X = xez[0];
				E = eta;
				Z = zeta;
			} else if (i == 1) {
				X = xi;
				E = xez[1];
				Z = zeta;
			} else if (i == 2) {
				X = xi;
				E = eta;
				Z = xez[2];
			} else if (i == 3) {
				X = xez[0];
				E = xez[1];
				Z = zeta;
			} else if (i == 4) {
				X = xi;
				E = xez[1];
				Z = xez[2];
			} else if (i == 5) {
				X = xez[0];
				E = eta;
				Z = xez[2];
			} else if (i == 6) {
				X = xez[0];
				E = xez[1];
				Z = xez[2];
			}
			tempVector.setX(0, X * h.component(0, 0) + E
					* h.component(1, 0) + Z
					* h.component(2, 0));
			tempVector.setX(1, X * h.component(0, 1) + E
					* h.component(1, 1) + Z
					* h.component(2, 1));
			tempVector.setX(2, X * h.component(0, 2) + E
					* h.component(1, 2) + Z
					* h.component(2, 2));

			refDistance = Math.sqrt(tempVector.squared());

			if (refDistance < oldRefDistance) {
				oldRefDistance = refDistance;
				nearestDr.E(tempVector);
			}
		} //end of for loop

		return nearestDr;
	}//end or compareDistance

	public void deform(etomica.space.Tensor deformationTensor) {
        for(int i=0; i<edgeVectors.length; i++) {
            edgeVectors[i].transform(deformationTensor);
        }
		update();
	}

	//        /**
	//         * Sets the length of each boundary edge to the corresponding value in the
	//         * given vector, keeping the shape of the box unchanged (other than the
	// change in size).
	//         */
	public void setDimensions(etomica.space.Vector v) {
		tempTensor.E(0);
		tempTensor.setComponent(0, 0, v.x(0));
		tempTensor.setComponent(1, 1, v.x(1));
		tempTensor.setComponent(2, 2, v.x(2));
		setDimensions(tempTensor);
	}

	/**
	 * Sets the boundary tensor to equal the given tensor.
	 */
	public void setDimensions(etomica.space.Tensor t) {
		h.E(t);
		update();
	}

	public double volume() {
		return Math.abs(h.component(0, 0)
				* h.component(1, 1)
				* h.component(2, 2));
	}

	/**
	 * imageOrigins and getOverFlowShifts are both probably incorrect, if they
	 * are even completed. They should definitely be checked before being
	 * implemented.
	 */

	int shellFormula, nImages, i, j, k, m;

	double[][] origins;

	public double[][] imageOrigins(int nShells) {
        int shellFormula = (2 * nShells) + 1;
        int nImages = space.powerD(shellFormula) - 1;
        double[][] origins = new double[nImages][space.D()];
        indexIterator.setSize(shellFormula);
        indexIterator.reset();
        int k = 0;
        while(indexIterator.hasNext()) {
            int[] index = indexIterator.next();
            temp1.E(index);
            temp1.PE(-(double)nShells);
            if(temp1.isZero()) continue;
            temp2.E(0.0);
            for(int i=0; i<space.D(); i++) {
                temp2.PEa1Tv1(temp1.x(i),edgeVectors[i]);
            }
            temp2.assignTo(origins[k++]);
        }
        return origins;
	}//end of imageOrigins

	//getOverflowShifts ends up being called by the display routines quite
	// often
	//so, in the interest of speed, i moved these outside of the function;

	public float[][] getOverflowShifts(etomica.space.Vector rr, double distance) {
		throw new RuntimeException(
				"Space3D.BoundaryDeformablePeriodic.getOverflowShifts not implmented");
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
    protected final boolean[] isPeriodic;
    private Tensor h;
    private Tensor hCopy;
    private final Vector[] edgeVectors; 
    private final Tensor hInv;
    private final Vector temp1, temp2;
    private final Vector dimensions;
    private final Vector dimensionsCopy;
    private final Vector dimensionsHalf;
    private final Vector tempVector;
    private final Vector nearestDr;
    private final Vector unit;
    private final Tensor tempTensor;
    private double oldRefDistance = 0.0;
    private double refDistance = 0.0;
    private final double[] workXEZ = new double[3];
    public boolean constantStrain = true;
    private boolean useBruteForceNearestImageMethod = true;
    private final IndexIteratorSequential indexIterator;

    // in the above, true uses brute force method!!

}//end of BoundaryDeformablePeriodic
