/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeActionTranslateTo;
import etomica.action.MoleculeChildAtomAction;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.space.Vector;
import etomica.atom.MoleculePair;
import etomica.space.Space;
import etomica.units.Degree;

/**
 *  Determine the second derivative of the atomic/ molecular potential energy w.r.t. to
 *   its generalized coordinates, u, where u is defined as the relative deviation of the 
 *   atom/ molecule from its nominal position 
 * 
 *  The class use <CoordinateDefinition> class SetToU method to put the atom/ molecule in
 *   space and the calculate the change in potential energy. 
 *   
 *  Output is the second derivative of the energy, phi: d2_phi/du_i^2, or d2_phi/(du_i.du_j)
 *  For Hessian Matrix:   d2_phi/(du_i.du_j) = d2_phi/(du_j.du_i)
 * 
 * 	This class is doing HALF ANALYTIC and HALF NUMERIC second derivative for self term
 * 	which, the first derivative is calculated using analytic approach and then compute the 
 *  second derivative using numerical derivative.
 * 	
 * @author taitan
 *
 */
public class CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen{
	
	public CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen(Space space, Box box, P2Nitrogen potential, CoordinateDefinitionNitrogen coordinateDefinition, boolean isAlpha){
		this(space, box, potential, coordinateDefinition, false, potential.getRange(), isAlpha);
	}
	
	public CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen(Space space, Box box, P2Nitrogen potential, CoordinateDefinitionNitrogen coordinateDefinition,
                                                            boolean doLatticeSum, double rC, boolean isAlpha){
		this.coordinateDefinition = coordinateDefinition;
		this.potential = potential;
		this.doLatticeSum = doLatticeSum;
		this.rC = rC;
		this.isAlpha = isAlpha;
		
		if(doLatticeSum){
			potential.setEnablePBC(false);
			potential.setRange(rC);
		}
		
		translateBy = new AtomActionTranslateBy(coordinateDefinition.getPrimitive().getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy); 
        pos = new MoleculePositionGeometricCenter(coordinateDefinition.getPrimitive().getSpace());
        translator = new MoleculeActionTranslateTo(coordinateDefinition.getPrimitive().getSpace());
        translator.setAtomPositionDefinition(pos);
        
		lsPosition = space.makeVector();
		destination = space.makeVector();
		workVec = space.makeVector();
        
        a = new double[ntab][ntab];
		generalizedCoord = new double[5];

		if(isAlpha){
			xVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(0);
			yVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(1);
			zVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(2); 
		} else {
			xVecBox = Math.sqrt(box.getBoundary().getEdgeVector(0).squared());
			yVecBox = Math.sqrt(box.getBoundary().getEdgeVector(1).squared());
			zVecBox = Math.sqrt(box.getBoundary().getEdgeVector(2).squared());
		}
	}
 	
	public double df(int[] moleculei, int dxi, double[] newU) {
		
		double[] u = newU;
		
		coordinateDefinition.setToUMoleculei(moleculei[0], newU);
		
		double rX = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(0);
		int nLayer = (int)(rC/rX + 0.5);
		MoleculePair pair = new MoleculePair();
		double df = 0.0;
		
		IMoleculeList moleculeList = coordinateDefinition.getBox().getMoleculeList();
		int numMolecule = moleculeList.getMoleculeCount();
		
		pair.atom0 = moleculeList.getMolecule(moleculei[0]);
		
		IMolecule molecule1;
		Vector[][] gradTorq;
		int firstRowMol = (int)Math.pow(numMolecule/1.99999, 1.0/3.0)*2;
		
		for (int i=0; i<numMolecule; i++){
			
			if(i==moleculei[0] && !doLatticeSum) continue;
			
			molecule1 = moleculeList.getMolecule(i); 
			pair.atom1 = molecule1;
			
			if(doLatticeSum){
				
				if(i==moleculei[0]){
					//THIS IS A HACK!
					int molNum;
					
					// alpha
					if(isAlpha){
						if((numMolecule-i)<5){
							molNum = (i-4);
						} else{
							molNum = (i+4);
						}
					} else {
						//beta
//						System.out.println("firstRowMol: " + firstRowMol);
						if(i%firstRowMol >=(firstRowMol-2)){
							molNum = (i-2);
						} else {
							molNum = (i+2);
						}
					}
					
					molecule1 = moleculeList.getMolecule(molNum);
					pair.atom1 = molecule1;

					//rotate the "borrowed" molecule 
					coordinateDefinition.setToUMoleculei(molNum, u);
					
					destination.E(pos.position(pair.atom0));
					translator.setDestination(destination);
					translator.actionPerformed(pair.atom1); 
				}
				
				if(isAlpha){
				
					for(int x=-nLayer; x<=nLayer; x++){
						for(int y=-nLayer; y<=nLayer; y++){
							for(int z=-nLayer; z<=nLayer; z++){
								if(i==moleculei[0] && x==0 && y==0 && z==0) continue;
	
								lsPosition.E(new double[]{x*xVecBox, y*yVecBox, z*zVecBox});
								translateBy.setTranslationVector(lsPosition);
								atomGroupActionTranslate.actionPerformed(molecule1);
	
								gradTorq = potential.gradientAndTorque(pair);
								
								if(dxi < 3){
									if(i!=moleculei[0]){
										df += gradTorq[0][0].getX(dxi);
									}
								
								} else{ // rotation dof
									workVec.E(gradTorq[1][0]);
									int rotAxisi = dxi==3 ? 2: 1;
									
									double dot = workVec.dot(coordinateDefinition.getMoleculeOrientation(pair.atom0)[rotAxisi]);
									if(dxi==3) dot*= -1;
									if(i==moleculei[0]) dot *= 2;
									df += dot;
	
								}
	
								lsPosition.TE(-1);
								translateBy.setTranslationVector(lsPosition);
								atomGroupActionTranslate.actionPerformed(molecule1);
							
							}	
						}	
					}
				} else {
					for(int x=-nLayer; x<=nLayer; x++){
						for(int y=-nLayer; y<=nLayer; y++){
							for(int z=-nLayer; z<=nLayer; z++){
								if(i==moleculei[0] && x==0 && y==0 && z==0) continue;
								lsPosition.E(new double[]{x*xVecBox-y*yVecBox*Math.cos(Degree.UNIT.toSim(60)), 
										 y*yVecBox*Math.sin(Degree.UNIT.toSim(60)), 
										 z*zVecBox});
								translateBy.setTranslationVector(lsPosition);
								atomGroupActionTranslate.actionPerformed(molecule1);
	
								gradTorq = potential.gradientAndTorque(pair);
								
								if(dxi < 3){
									if(i!=moleculei[0]){
										df += gradTorq[0][0].getX(dxi);
									}
								
								} else{ // rotation dof
									workVec.E(gradTorq[1][0]);
									int rotAxisi = dxi==3 ? 2: 1;
									
									double dot = workVec.dot(coordinateDefinition.getMoleculeOrientation(pair.atom0)[rotAxisi]);
									if(dxi==3) dot*= -1;
									if(i==moleculei[0]) dot *= 2;
									df += dot;
	
								}
	
								lsPosition.TE(-1);
								translateBy.setTranslationVector(lsPosition);
								atomGroupActionTranslate.actionPerformed(molecule1);
							
							}	
						}	
					}
					
				}
				
				//PUTTING the HACK molecule back to its initial position
				if(i==moleculei[0]){
					int molNum;
					
					//alpha
					if(isAlpha){
						if((numMolecule-i)<5){
							molNum = (i-4);
						} else{
							molNum = (i+4);
						}
					} else {
						//beta
						
						if(i%firstRowMol >=(firstRowMol-2)){
							molNum = (i-2);
						} else {
							molNum = (i+2);
						}
					}
					for(int y=0; y<u.length; y++){
						u[y] = 0.0;
					}
					// putting the molecule back its own lattice site
					// and rotate back to original orientation
					coordinateDefinition.setToUMoleculei(molNum, u);
				}
				
			} else 	{   
				gradTorq = potential.gradientAndTorque(pair);
				
				// Cartesian x-, y- and z-
				if(dxi < 3){
					df += gradTorq[0][0].getX(dxi);
				
				} else{ // rotation dof
					workVec.E(gradTorq[1][0]);
					int rotAxisi = dxi==3 ? 2: 1;
					
					double dot = workVec.dot(coordinateDefinition.getMoleculeOrientation(pair.atom0)[rotAxisi]);
					if(dxi==3) dot*= -1;
					df += dot;	
				}
			}
			
		}
		return df;
		
	}
	
	public double d2phi_du2(int[] moleculei, int[] d) {
		
		
		if(moleculei[0] != moleculei[1]){
			throw new RuntimeException("<CalcHalfAnalyticHalfNumeric2ndDerivativeNitrogen>: "
					+ " only do self-term!! for non-selfterm, please refer to <CalcAnalytical2ndDerivativeNitrogen>");
		}
		/*
		 * d[i] = the i-th elements of the derivative
		 * u[i] = the generalized coordinate/ relative displacement from 
		 * 			nominal position 
		 */
        deltaU = fixedDeltaU;
        a[0][0] = computeA(moleculei, d, deltaU);
        
		double err = big;
        double d2fdu2 = Double.NaN;

		for(int i=1; i<ntab; i++){
			deltaU = deltaU /con;
			
            a[0][i] = computeA(moleculei, d, deltaU);
			
            fac = con2;
			
			for(int j=1; j<=i; j++){
				a[j][i] = (a[j-1][i]*fac - a[j-1][i-1])/(fac-1);
				fac = con2*fac;
				errt = Math.max(Math.abs(a[j][i]-a[j-1][i]), Math.abs(a[j][i]-a[j-1][i-1]));
				
				if (errt <= err){
					err = errt;
					d2fdu2 = a[j][i];
				}
			}
			
			if (Math.abs(a[i][i]-a[i-1][i-1]) >= safe*err){
				break;
			}
		}

		return d2fdu2;
	}

	public double computeA(int[] moleculei, int[] d, double deltaU){
		
		if(d.length > 2){
			throw new RuntimeException("In computeA method; Derivative method can only do first and" +
					"second order derivative!!!");
		}
	
		/*
		 *  Doing FIRST DERIVATIVE of a First Derivative
		 *  
		 *  f is the energy function
		 */
		
			/*
			 * 
			 *  when u_i = u_j
			 *  
			 *      df                f(u_i + deltaU) - f(u_i - deltaU) 
			 *  -----------  =   ----------------------------------------------- 
			 *     du_i                           2*deltaU
			 *  
			 */
				
			generalizedCoord[d[1]] = deltaU;
			double f_up1 = df(moleculei, d[0], generalizedCoord);
	
			generalizedCoord[d[1]] = -deltaU;
			double f_um1 = df(moleculei, d[0], generalizedCoord);
				
			generalizedCoord[d[1]] = 0.0;
				
			setToInitialPosition(moleculei);
			return (f_up1 - f_um1)/(2*deltaU);  
		
	}
	
	public void setToInitialPosition(int[] moleculei){
		double newU[] = new double[5];
		coordinateDefinition.setToUMoleculei(moleculei[0], newU);
		coordinateDefinition.setToUMoleculei(moleculei[1], newU);

	}
	
	public double getFixedDeltaU() {
		return fixedDeltaU;
	}

	public void setFixedDeltaU(double fixedDeltaU) {
		this.fixedDeltaU = fixedDeltaU;
	}

	protected Box box;
	protected MoleculePositionGeometricCenter pos;
	protected MoleculeActionTranslateTo translator;
	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected P2Nitrogen potential;
	protected AtomActionTranslateBy translateBy;
	protected MoleculeChildAtomAction atomGroupActionTranslate;
	protected Vector lsPosition, destination, workVec;
	protected double errt, fac, xVecBox, yVecBox, zVecBox, rC;
	protected double deltaU;
	protected double [][] a;
	protected double [] generalizedCoord;
	protected boolean doLatticeSum = false;
	double fixedDeltaU = 0.01;
	final int ntab = 10;
	final double con = 1.4;
	final double con2 = con*con;
	final double big = Double.MAX_VALUE;
	final double safe = 2.0;
	final boolean isAlpha;
	
}
