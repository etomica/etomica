package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.MoleculePair;
import etomica.data.types.DataTensor;
import etomica.space.ISpace;
import etomica.space.Space;

/**
 *  Determine the second derivative of the atomic/ molecular potential energy w.r.t. to
 *   its generalized coordinates, u, where u is defined as the relative deviation of the 
 *   atom/ molecule from its nominal position 
 * 
 *  The class use <CoordinateDefinition> class SetToU method to put the atom/ molecule in
 *   space and the calculate the change in potential energy. 
 *   
 * 
 * @author taitan
 *
 */
public class CalcAnalytical2ndDerivativeNitrogen{
	
	public CalcAnalytical2ndDerivativeNitrogen(ISpace space, IBox box, P2Nitrogen potential,CoordinateDefinitionNitrogen coordinateDefinition){
		this(space, box, potential, coordinateDefinition, false, potential.getRange());
	}
	
	public CalcAnalytical2ndDerivativeNitrogen(ISpace space, IBox box, P2Nitrogen potential,CoordinateDefinitionNitrogen coordinateDefinition,
			boolean doLatticeSum, double rC){
		this.coordinateDefinition = coordinateDefinition;
		this.potential = potential;
		this.doLatticeSum = doLatticeSum;
		this.rC = rC;
		this.space = space;
		
		this.potential = potential;
		if(doLatticeSum){
			potential.setEnablePBC(false);
			potential.setRange(rC);
		}
		
		translateBy = new AtomActionTranslateBy(coordinateDefinition.getPrimitive().getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy); 
        pos = new AtomPositionGeometricCenter(coordinateDefinition.getPrimitive().getSpace());
        translator = new MoleculeActionTranslateTo(coordinateDefinition.getPrimitive().getSpace());
        translator.setAtomPositionDefinition(pos);
        
		lsPosition = space.makeVector();
		destination = space.makeVector();
		com1 = space.makeVector();
		com2 = space.makeVector();
		dr = space.makeVector();
		
		xVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(0);
		yVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(1);
		zVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(2); 
		
		secDerXr = new IVectorMutable[2][3];
		for(int i=0; i<secDerXr.length; i++){
			for(int j=0; j<secDerXr[0].length; j++){
				secDerXr[i][j] = space.makeVector(); 
			}	
		}
		
		dUdRotA = new IVectorMutable[2];
		dUdRotB = new IVectorMutable[2];
		
		for(int i=0; i<dUdRotA.length; i++){
			dUdRotA[i] = space.makeVector();
			dUdRotB[i] = space.makeVector(); 
		}	
		
		int numMolec = coordinateDefinition.getBox().getMoleculeList().getMoleculeCount();
		initMolecOrientation = new IVectorMutable[numMolec][3];
		
		for (int i=0; i<numMolec; i++){
			initMolecOrientation[i] = space.makeVectorArray(3);
			initMolecOrientation[i] = coordinateDefinition.getMoleculeOrientation(box.getMoleculeList().getMolecule(i));
		}
	}
 	
	public double[][] d2phi_du2(int[] moleculei) {
		
		if(moleculei[0] == moleculei[1]) {
			throw new RuntimeException("C<CalcAnalytical2ndDerivcationNitrogen> CANNOT HANDLE SELF-TERM YET!");
		}
		
		double [][] d2r = new double[5][5];
		
		MoleculePair pair = new MoleculePair();
		pair.atom0 = coordinateDefinition.getBox().getMoleculeList().getMolecule(moleculei[0]);
		pair.atom1 = coordinateDefinition.getBox().getMoleculeList().getMolecule(moleculei[1]);
		
		DataTensor tensorTrans = new DataTensor(Space.getInstance(3));
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
//		IMolecule nitrogena = pair.getMolecule(0);
//		IMolecule nitrogenb = pair.getMolecule(1);
//		
//		IVectorMutable pos1 = (nitrogena.getChildList().getAtom(1)).getPosition();
//		IVectorMutable pos2 = (nitrogenb.getChildList().getAtom(1)).getPosition();
//		
//		com1.E(pos1);
//		com2.E(pos2);
//		
//		IVectorMutable diff1 = space.makeVector();
//		IVectorMutable diff2 = space.makeVector();
//		
//		diff1.Ev1Mv2(com1, nitrogena.getChildList().getAtom(0).getPosition());
//		diff2.Ev1Mv2(com2, nitrogenb.getChildList().getAtom(0).getPosition());
//					
//		com1.PEa1Tv1(-0.5, diff1); 		
//		com2.PEa1Tv1(-0.5, diff2);
//		
//		// Vertical 3x2 matrix  (upper right corner)
//		//dudrU3MolA
//		dUdRotA[0].E(new double[]{ d2r[0][3], d2r[1][3], d2r[2][3] }); 
//		//dudrU4MolA
//		dUdRotA[1].E(new double[]{ d2r[0][4], d2r[1][4], d2r[2][4] }); 
//		
//		
//		// Horizontal 2x3 matrix (bottom left corner)
//		//dudrU3MolB
//		dUdRotB[0].E(new double[]{ d2r[3][0], d2r[3][1], d2r[3][2] }); 
//		//dudrU4MolB
//		dUdRotB[1].E(new double[]{ d2r[4][0], d2r[4][1], d2r[4][2] }); 
//		
//		IVectorMutable workVec = space.makeVector();
//		IVectorMutable sumVec = space.makeVector();
//		
//		for(int i=0; i<2; i++){
//			sumVec.E(0.0);
//			for(int iAtomNB=0; iAtomNB<nitrogenb.getChildList().getAtomCount(); iAtomNB++){
//				dr.Ev1Mv2(nitrogenb.getChildList().getAtom(iAtomNB).getPosition(), com2);
//				workVec.E(dUdRotA[i]);
//				workVec.XE(dr);
//				sumVec.PE(workVec);
//			}
//			d2r[3][3+i] = sumVec.dot(initMolecOrientation[][]);
//		}
		
		return d2r;
	}
	
	protected IVectorMutable[][] initMolecOrientation;
	protected IVectorMutable[][] secDerXr;
	protected IVectorMutable[] dUdRotA, dUdRotB;
	protected IBox box;
	protected ISpace space;
	protected AtomPositionGeometricCenter pos;
	protected MoleculeActionTranslateTo translator;
	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected P2Nitrogen potential;
	protected AtomActionTranslateBy translateBy;
	protected MoleculeChildAtomAction atomGroupActionTranslate;
	protected IVectorMutable lsPosition, destination, com1, com2, dr;
	protected double xVecBox, yVecBox, zVecBox, rC;
	protected boolean doLatticeSum = false;
	
	
	

	public void test(IMoleculeList pair){
	

	
	    
	  						        

	}
	
	
}
