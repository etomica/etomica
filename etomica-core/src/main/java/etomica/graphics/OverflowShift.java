/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.space.BoundaryRectangular;
import etomica.space.Space;
import etomica.space3d.BoundaryTruncatedOctahedron;

public class OverflowShift {

    private int dim;

	OverflowShift(Space space) {
		dim = space.D();

		
	}
	
	/**
	 * Provides information needed so that drawing can be done of portions of
	 * the periodic images of an atom that overflow into the volume (because the
	 * periodic image is just outside the volume).
	 * 
	 * @param boundary
	 *            boundary of atoms
	 * @param rr
	 *            position of the atom in the central image
	 * @param distance
	 *            size of the atom, indicating how far its image extends
	 * @return all displacements needed to draw all overflow images; first index
	 *         indicates each displacement, second index is the xyz translation
	 *         needed to the overflow image
	 */
    public float[][] getShifts(Boundary boundary, Vector rr, double distance) {
    	
		boolean needShift[] = new boolean[dim];
		float[][] shifts = null;
		Vector dimensions = boundary.getBoxSize();

		//
		// BoundaryRectangular
		//
		if(boundary instanceof BoundaryRectangular) {

	        int numVectors = 1;
	        for (int i=0; i<dim; i++) {
	           if ((rr.getX(i) - distance < -0.5*dimensions.getX(i)) || (rr.getX(i) + distance > 0.5*dimensions.getX(i))) {
	              //each previous vector will need an additional copy in this dimension 
	              numVectors *= 2;
	              //remember that
	              needShift[i] = true;
	           }
	           else {
	              needShift[i] = false;
	           }
	        }
	        
	        if(numVectors == 1) {
	        	shifts = new float[0][0];
	        }
	        else {
		        shifts = new float[numVectors][dim];
		        double[] rrArray = new double[dim];
		        rr.assignTo(rrArray);
		        for (int i=0; i<dim; i++) {
		//           shifts[0][i] = (float)rrArray[i];
		        }
		        int iVector = 1;
		
		        for (int i=0; iVector<numVectors; i++) {
		           if (!needShift[i]) {
		              //no shift needed for this dimension
		              continue;
		           }
		           double delta = -dimensions.getX(i);
		           if (rr.getX(i) - distance < -0.5*dimensions.getX(i)) {
		              delta = -delta;
		           }
		           //copy all previous vectors and apply a shift of delta to the copies
		           for (int jVector=0; jVector<iVector; jVector++) {
		              for (int j=0; j<dim; j++) {
		                  shifts[jVector+iVector][j] = shifts[jVector][j];
		              }
		              shifts[jVector+iVector][i] += delta;
		           }
		           iVector *= 2;
		        }
	        }
		}
		//
		// BoundaryTruncatedOctahedron
		//
		else if(boundary instanceof BoundaryTruncatedOctahedron) {
	        int numVectors = 1;
	        for (int i = 0; i < dim; i++) {
	            if ((rr.getX(i) - distance < 0.0)
	                    || (rr.getX(i) + distance > dimensions.getX(i))) {
	                //each previous vector will need an additional copy in this
	                // dimension
	                numVectors *= 2;
	                //remember that
	                needShift[i] = true;
	            } else {
	                needShift[i] = false;
	            }
	        }

	        if (numVectors == 1) {
	        	shifts = new float[0][0];
	        }
	        else {
		        shifts = new float[numVectors][dim];
		        double[] rrArray = new double[3];
		        rr.assignTo(rrArray);
		        for (int i = 0; i < dim; i++) {
		            shifts[0][i] = (float) rrArray[i];
		        }
		        int iVector = 1;
	
		        for (int i = 0; iVector < numVectors; i++) {
		            if (!needShift[i]) {
		                //no shift needed for this dimension
		                continue;
		            }
		            double delta = -dimensions.getX(i);
		            if (rr.getX(i) - distance < 0.0) {
		                delta = -delta;
		            }
		            //copy all previous vectors and apply a shift of delta to the
		            // copies
		            for (int jVector = 0; jVector < iVector; jVector++) {
		                for (int j = 0; j < dim; j++) {
		                    shifts[jVector + iVector][j] = shifts[jVector][j];
		                }
		                shifts[jVector + iVector][i] += delta;
		            }
		            iVector *= 2;
		        }
	        }
		}
		else {
			shifts = new float[0][0];
		}

        return shifts;
     }
     
}
