/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.box.Box;
import etomica.data.types.DataTensor;
import etomica.molecule.MoleculePair;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;

/**
 *  Determine the second derivative of the atomic/ molecular potential energy w.r.t. to
 *   its generalized coordinates, u, where u is defined as the relative deviation of the 
 *   atom/ molecule from its nominal position 
 * 	
 *  Analytic expression of the derivative for interparticle interaction ONLY
 *  - this DOES NOT include the self-term
 * 
 * @author Tai Boon Tan & Andrew Schultz
 *
 */
public class CalcAnalytical2ndDerivativeNitrogen{
	
	public CalcAnalytical2ndDerivativeNitrogen(Space space, Box box, P2Nitrogen potential, CoordinateDefinitionNitrogen coordinateDefinition){
		this(space, box, potential, coordinateDefinition, false, potential.getRange());
	}
	
	public CalcAnalytical2ndDerivativeNitrogen(Space space, Box box, P2Nitrogen potential, CoordinateDefinitionNitrogen coordinateDefinition,
                                               boolean doLatticeSum, double rC){
		this.coordinateDefinition = coordinateDefinition;
		this.potential = potential;
		this.doLatticeSum = doLatticeSum;
		this.space = space;
		
		this.potential = potential;
		if(doLatticeSum){
			potential.setEnablePBC(false);
			potential.setRange(rC);
		}
		
		workVec = space.makeVector();
		
		secDerXr = new Vector[2][3];
		for(int i=0; i<secDerXr.length; i++){
			for(int j=0; j<secDerXr[0].length; j++){
				secDerXr[i][j] = space.makeVector(); 
			}	
		}
		
		dUdRotA = new Vector[2];
		dUdRotB = new Vector[2];
		
		for(int i=0; i<dUdRotA.length; i++){
			dUdRotA[i] = space.makeVector();
			dUdRotB[i] = space.makeVector(); 
		}	
		
		int numMolec = coordinateDefinition.getBox().getMoleculeList().getMoleculeCount();
		initMolecOrientation = new Vector[numMolec][3];
		
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDefinition.getMoleculeOrientation(box.getMoleculeList().getMolecule(i));
		}
	}
 	
	public double[][] d2phi_du2(int[] moleculei) {
		
		if(moleculei[0] == moleculei[1]) {
			throw new RuntimeException("C<CalcAnalytical2ndDerivcationNitrogen> CANNOT HANDLE SELF-TERM YET!");
		}
		
		pair.atom0 = coordinateDefinition.getBox().getMoleculeList().getMolecule(moleculei[0]);
		pair.atom1 = coordinateDefinition.getBox().getMoleculeList().getMolecule(moleculei[1]);
		
		tensorTrans.E(potential.secondDerivative(pair));
		
		// i-row and j-column
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				d2r[i][j] = tensorTrans.x.component(i, j);
			}	
		}
		
		secDerXr = potential.secondDerivativeXr(pair);
		
		//Filling up 4th and 5th column (Rotation moleculeA) 
		//       (1st, 2nd and 3rd Row) (Translation molecule B)
		// when j=3, rotation in u3; j=4, rotation in u4
		// the negative is to account for the opposite direction (direction of the torque)	
		for(int i=0; i<3; i++){
			d2r[i][3] = -initMolecOrientation[moleculei[0]][2].dot(secDerXr[0][i]);
			d2r[i][4] = initMolecOrientation[moleculei[0]][1].dot(secDerXr[0][i]);
			
		}	
		
		//Filling up 4th and 5th row (Rotation moleculeB) 
		//       (1st, 2nd and 3rd column) (Translation molecule A)
		// when j=3, rotation in u3; j=4, rotation in u4
		for(int i=0; i<3; i++){
			d2r[3][i] = -initMolecOrientation[moleculei[1]][2].dot(secDerXr[1][i]); 
			d2r[4][i] = initMolecOrientation[moleculei[1]][1].dot(secDerXr[1][i]); 
			
		}	
				
		Tensor tensor = potential.secondDerivativeXrRotRot(pair);
	
		for(int i=0; i<2; i++){
			for(int j=0; j<2; j++){
					
				workVec.E(initMolecOrientation[moleculei[1]][2-j]);
				tensor.transform(workVec);

				d2r[3+j][3+i] = workVec.dot(initMolecOrientation[moleculei[0]][2-i]);
				if(i!=j) d2r[3+j][3+i] *= -1;
				
			}
		}
		
		return d2r;
	}
	
	protected Vector[][] initMolecOrientation;
	protected Vector[][] secDerXr;
	protected Vector[] dUdRotA, dUdRotB;
	protected Box box;
	protected Space space;
	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected P2Nitrogen potential;
	protected Vector workVec;
	protected boolean doLatticeSum = false;
    protected double [][] d2r = new double[5][5];
    protected MoleculePair pair = new MoleculePair();
    protected DataTensor tensorTrans = new DataTensor(Space.getInstance(3));
}
