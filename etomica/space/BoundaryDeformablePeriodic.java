package etomica.space;

import etomica.math.geometry.Parallelepiped;
import etomica.math.geometry.Polytope;
import etomica.simulation.Simulation;
import etomica.space3d.Vector3D;

/**
 * @author skkwak
 */

//nan needs a cleanup of imageOrigins & getOverflowShifts.
public class BoundaryDeformablePeriodic extends Boundary implements BoundaryPeriodic {

	protected final boolean[] isPeriodic;
	private Tensor boundaryTensor;
	private Tensor boundaryTensorCopy;
	private final Vector temp;
	private final Vector dimensions;
	private final Vector dimensionsCopy;
	private final Vector dimensionsHalf;
	private final Vector workVector;
	private final Vector nearestDr;
	private final Tensor workTensor;
	private final Tensor tempT;
	private double oldRefDistance = 0.0;
	private double refDistance = 0.0;
	private final double[] workXEZ = new double[3];
	public boolean constantStrain = true;
	private boolean useBruteForceNearestImageMethod = true; 
	// in the above, true uses brute force method!!

	public BoundaryDeformablePeriodic(Simulation sim) {
		this(sim.space, makePeriodicity(sim.space.D()), sim.getDefaults().boxSize);
	}

	public BoundaryDeformablePeriodic(Space space, boolean[] periodicity, double boxSize) {
	    //nan 2D
	    this(space, periodicity, boxSize, new Vector3D (1,0,0), new Vector3D (0,1,0), 
	            new Vector3D (0,0,1));	
	}
	
//	nan 2D
	public BoundaryDeformablePeriodic(Space space, boolean[] periodicity, double boxSize, 
	        Vector3D a, Vector3D b, Vector3D c) {
	    //super(space, new Parallelepiped(space, a, b, c));
//	  nan 2D
		super(space, makeShape(space, a, b, c));
		isPeriodic = (boolean[]) periodicity.clone();
		boundaryTensor = space.makeTensor();
		boundaryTensorCopy = space.makeTensor();
		temp = space.makeVector();
		dimensions = space.makeVector();
		dimensionsCopy = space.makeVector();
		dimensionsHalf = space.makeVector();
		workVector = space.makeVector();
		nearestDr = space.makeVector();
		workTensor = space.makeTensor();
		tempT = space.makeTensor();
		workVector.E(boxSize);
		setDimensions(workVector);
		makeStrainTensor();
	}
//	nan 2D
    private static Polytope makeShape(Space space, Vector a, Vector b, Vector c) {
        switch(space.D()) {
            case 1: throw new IllegalArgumentException("Wrong number of dimensions for BoundaryDeformablePeriodic");
            case 2: throw new IllegalArgumentException("Wrong number of dimensions for BoundaryDeformablePeriodic");
            case 3: return new Parallelepiped(space, (Vector3D)a, (Vector3D)b,
                    (Vector3D)c);
            default: throw new IllegalArgumentException("BoundaryRectangular not appropriate to given space");
        }
    }
	public boolean[] getPeriodicity() {
		return isPeriodic;
	}

	private static boolean[] makePeriodicity(int D) {
		boolean[] isPeriodic = new boolean[D];
		for (int i = 0; i < D; i++) {
			isPeriodic[i] = true;
		}
		return isPeriodic;
	}

	public final etomica.space.Vector dimensions() {
		return dimensionsCopy;
	}

	public etomica.space.Tensor boundaryTensor() {
		return boundaryTensorCopy;
	}

	public etomica.space.Vector randomPosition() {
		temp.setRandom(1.0);
		temp.transform(boundaryTensor);
		return temp;
	}

	private void updateCopies() {
		dimensionsHalf.Ea1Tv1(0.5, dimensions);
		dimensionsCopy.E(dimensions);
		boundaryTensorCopy.E(boundaryTensor);
	}

	public void nearestImage(Vector dr) {
		checkNearestImage(dr);
	}

	public Vector centralImage(Vector r) {
		temp.E(r);
		checkCentralImage(temp);
		temp.ME(r);
		return temp;
	}

	public void setConstantStrain(boolean b) {
		constantStrain = b;
	}

	/**
	 * @return Returns the constantStrain.
	 */
	public boolean isConstantStrain() {
		return constantStrain;
	}
	
	public void makeStrainTensor() {
		tempT.E(boundaryTensorCopy);
		tempT.transpose(); //r=h*s -->h matrix is a transpose of boundaryTensor
		tempT.inverse(); //to get s, inverse of h times r, which is a position
						 // of atom
	}

	public void setUseBruteForceNearestImageMethod(boolean bb) {
		useBruteForceNearestImageMethod = bb;
	}
//skkwak next!!!
	public void checkNearestImage(Vector dr2) { // dr2 is not dr^2
		if (!useBruteForceNearestImageMethod) {
			if (!constantStrain) {
				makeStrainTensor();
			} //in Main Class, makeStrainTensor() must be called then set up
			  // boolean of contantStrain!!!

			double xi = dr2.x(0) * tempT.component(0, 0) + dr2.x(1)
					* tempT.component(0, 1) + dr2.x(2) * tempT.component(0, 2);
			double eta = dr2.x(0) * tempT.component(1, 0) + dr2.x(1)
					* tempT.component(1, 1) + dr2.x(2) * tempT.component(1, 2);
			double zeta = dr2.x(0) * tempT.component(2, 0) + dr2.x(1)
					* tempT.component(2, 1) + dr2.x(2) * tempT.component(2, 2);

			if (xi > -1.0e-10 && xi < 1.0e-10) {
				xi = 0;
			}
			if (eta > -1.0e-10 && eta < 1.0e-10) {
				eta = 0;
			}
			if (zeta > -1.0e-10 && zeta < 1.0e-10) {
				zeta = 0;
			}

			xi += 0.5;
			eta += 0.5;
			zeta += 0.5;

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

			xi -= 0.5;
			eta -= 0.5;
			zeta -= 0.5;

			dr2.setX(0, xi * boundaryTensor.component(0, 0) + eta
					* boundaryTensor.component(1, 0) + zeta
					* boundaryTensor.component(2, 0));
			dr2.setX(1, xi * boundaryTensor.component(0, 1) + eta
					* boundaryTensor.component(1, 1) + zeta
					* boundaryTensor.component(2, 1));
			dr2.setX(2, xi * boundaryTensor.component(0, 2) + eta
					* boundaryTensor.component(1, 2) + zeta
					* boundaryTensor.component(2, 2));

			//                    } //inner if
		} else {
			if (!constantStrain) {
				makeStrainTensor();
			} //in Main Class, makeStrainTensor() must be called then set up
			  // boolean of contantStrain!!!

			oldRefDistance = Math.sqrt(dr2.squared());
			nearestDr.E(dr2);

			double xi = dr2.x(0) * tempT.component(0, 0) + dr2.x(1)
					* tempT.component(0, 1) + dr2.x(2) * tempT.component(0, 2);
			double eta = dr2.x(0) * tempT.component(1, 0) + dr2.x(1)
					* tempT.component(1, 1) + dr2.x(2) * tempT.component(1, 2);
			double zeta = dr2.x(0) * tempT.component(2, 0) + dr2.x(1)
					* tempT.component(2, 1) + dr2.x(2) * tempT.component(2, 2);

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
			workVector.setX(0, X * boundaryTensor.component(0, 0) + E
					* boundaryTensor.component(1, 0) + Z
					* boundaryTensor.component(2, 0));
			workVector.setX(1, X * boundaryTensor.component(0, 1) + E
					* boundaryTensor.component(1, 1) + Z
					* boundaryTensor.component(2, 1));
			workVector.setX(2, X * boundaryTensor.component(0, 2) + E
					* boundaryTensor.component(1, 2) + Z
					* boundaryTensor.component(2, 2));

			refDistance = Math.sqrt(workVector.squared());

			if (refDistance < oldRefDistance) {
				oldRefDistance = refDistance;
				nearestDr.E(workVector);
			}
		} //end of for loop

		return nearestDr;
	}//end or compareDistance

	public void checkCentralImage(Vector r2) {
		if (!constantStrain) {
			makeStrainTensor();
		} //in Main Class, makeStrainTensor() must be called then set up
		  // boolean of contantStrain!!!
		double xi = r2.x(0) * tempT.component(0, 0) + r2.x(1)
				* tempT.component(0, 1) + r2.x(2) * tempT.component(0, 2);
		double eta = r2.x(0) * tempT.component(1, 0) + r2.x(1)
				* tempT.component(1, 1) + r2.x(2) * tempT.component(1, 2);
		double zeta = r2.x(0) * tempT.component(2, 0) + r2.x(1)
				* tempT.component(2, 1) + r2.x(2) * tempT.component(2, 2);

		if (Math.abs(xi) < 1e-10) {
			xi = 0;
		}
		if (Math.abs(eta) < 1e-10) {
			eta = 0;
		}
		if (Math.abs(zeta) < 1e-10) {
			zeta = 0;
		}

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

		r2.setX(0, xi * boundaryTensorCopy.component(0, 0) + eta
				* boundaryTensorCopy.component(1, 0) + zeta
				* boundaryTensorCopy.component(2, 0));
		r2.setX(1, xi * boundaryTensorCopy.component(0, 1) + eta
				* boundaryTensorCopy.component(1, 1) + zeta
				* boundaryTensorCopy.component(2, 1));
		r2.setX(2, xi * boundaryTensorCopy.component(0, 2) + eta
				* boundaryTensorCopy.component(1, 2) + zeta
				* boundaryTensorCopy.component(2, 2));

	}//end checkCentralImage

	private void updateDimensions() {
		for (i = 0; i < space.D(); i++)
			dimensions.setX(i, Math.abs(boundaryTensor.component(i, i)));
	}

	public void deform(etomica.space.Tensor deformationTensor) {
		boundaryTensor.TE(deformationTensor);
		updateDimensions();
		updateCopies();
		makeStrainTensor();
	}

	//        /**
	//         * Sets the length of each boundary edge to the corresponding value in the
	//         * given vector, keeping the shape of the box unchanged (other than the
	// change in size).
	//         */
	public void setDimensions(etomica.space.Vector v) {
		workTensor.E(0);
		workTensor.setComponent(0, 0, v.x(0));
		workTensor.setComponent(1, 1, v.x(1));
		workTensor.setComponent(2, 2, v.x(2));
		setDimensions(workTensor);
	}

	/**
	 * Sets the boundary tensor to equal the given tensor.
	 */
	public void setDimensions(etomica.space.Tensor t) {
		boundaryTensor.E(t);
		updateDimensions();
		updateCopies();
	}

	public double volume() {
		return Math.abs(boundaryTensor.component(0, 0)
				* boundaryTensor.component(1, 1)
				* boundaryTensor.component(2, 2));
	}

	/**
	 * imageOrigins and getOverFlowShifts are both probably incorrect, if they
	 * are even completed. They should definitely be checked before being
	 * implemented.
	 */

	int shellFormula, nImages, i, j, k, m;

	double[][] origins;

	public double[][] imageOrigins(int nShells) {
		throw new RuntimeException(
				"imageOrigins not implemented in Space3D.BoundaryDeformablePeriodic");
		/*
		 * shellFormula = (2 * nShells) + 1; nImages =
		 * shellFormula*shellFormula*shellFormula-1; origins = new
		 * double[nImages][D]; for (k=0,i=-nShells; i <=nShells; i++) { for
		 * (j=-nShells; j <=nShells; j++) { for (m=-nShells; m <=nShells; m++) {
		 * if ((i==0 && j==0) && m==0 ) {continue;} origins[k][0] =
		 * i*dimensions.x; origins[k][1] = j*dimensions.y; origins[k][2] =
		 * m*dimensions.z; k++; } } } return origins;
		 */
	}//end of imageOrigins

	//getOverflowShifts ends up being called by the display routines quite
	// often
	//so, in the interest of speed, i moved these outside of the function;
	int shiftX, shiftY, shiftZ;

	Vector r;

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

}//end of BoundaryDeformablePeriodic
