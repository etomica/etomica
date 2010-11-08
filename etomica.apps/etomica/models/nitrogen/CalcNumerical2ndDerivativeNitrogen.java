package etomica.models.nitrogen;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVectorMutable;
import etomica.atom.MoleculePair;

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
 *  NOTE: fixedDeltaU in the field is the interval of change. It has to be CAREFULLY picked!
 *        Usually, start from a larger value (~0.1, e.g.) and the algorithm will reduce 
 *        the value. One way to check is to compare d2f/(duiduj) and d2f/(dujdui), of course
 *        the smaller the value the better it is.
 * 
 * @author taitan
 *
 */
public class CalcNumerical2ndDerivativeNitrogen{
	
	public CalcNumerical2ndDerivativeNitrogen(IBox box, P2Nitrogen potential,CoordinateDefinitionNitrogen coordinateDefinition){
		this(box, potential, coordinateDefinition, false);
	}
	
	public CalcNumerical2ndDerivativeNitrogen(IBox box, P2Nitrogen potential,CoordinateDefinitionNitrogen coordinateDefinition,
			boolean doLatticeSum){
		this.coordinateDefinition = coordinateDefinition;
		this.potential = potential;
		this.doLatticeSum = doLatticeSum;
		
		if(doLatticeSum){
			potential.setEnablePBC(false);
		}
		translateBy = new AtomActionTranslateBy(coordinateDefinition.getPrimitive().getSpace());
        atomGroupActionTranslate = new MoleculeChildAtomAction(translateBy); 
		lsPosition = coordinateDefinition.getPrimitive().getSpace().makeVector();
        
        a = new double[ntab][ntab];
		generalizedCoord = new double[2][5];

		xVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(0);
		yVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(1);
		zVecBox = coordinateDefinition.getBox().getBoundary().getBoxSize().getX(2); 
		
	}
 	
	public double f(int[] moleculei, double[][] newU) {
		
		if(moleculei[0] == moleculei[1]){
			double[] u = new double[5];
			for(int i=0; i<newU[0].length; i++){
				u[i] += newU[0][i];
				u[i] += newU[1][i];
			}
			coordinateDefinition.setToUMoleculei(moleculei[0], u);
			
		} else{
			coordinateDefinition.setToUMoleculei(moleculei[0], newU[0]);
			coordinateDefinition.setToUMoleculei(moleculei[1], newU[1]);
		}
		
		MoleculePair pair = new MoleculePair();
		int nLayer= 10;
		double sum = 0.0;
		
		//	pair of identical molecules 
		if(moleculei[0] == moleculei[1]){
			IMoleculeList moleculeList = coordinateDefinition.getBox().getMoleculeList();
			int numMolecule = moleculeList.getMoleculeCount();
			
			pair.atom0 = moleculeList.getMolecule(moleculei[0]);
			for (int i=0; i<numMolecule; i++){
				if(i==moleculei[0]) continue; 
				IMolecule molecule1 = moleculeList.getMolecule(i); 
				pair.atom1 = molecule1;
				
				if(doLatticeSum){
					for(double x=-xVecBox*nLayer; x<=xVecBox*nLayer; x+=xVecBox){
						for(double y=-yVecBox*nLayer; y<=yVecBox*nLayer; y+=yVecBox){
							for(double z=-zVecBox*nLayer; z<=zVecBox*nLayer; z+=zVecBox){
								lsPosition.E(new double[]{x, y, z});
								translateBy.setTranslationVector(lsPosition);
								atomGroupActionTranslate.actionPerformed(molecule1);
			
								sum += potential.energy(pair);
								
								lsPosition.TE(-1);
								translateBy.setTranslationVector(lsPosition);
								atomGroupActionTranslate.actionPerformed(molecule1);
							}	
						}	
					}
					
				} else {
		
					sum += potential.energy(pair);
				}
				
			}
			return sum;
		}
		
		//pair of non identical molecules 
		pair.atom0 = coordinateDefinition.getBox().getMoleculeList().getMolecule(moleculei[0]);
		IMolecule molecule1 = coordinateDefinition.getBox().getMoleculeList().getMolecule(moleculei[1]); 
		pair.atom1 = molecule1;

		if(doLatticeSum){
			for(double x=-xVecBox*nLayer; x<=xVecBox*nLayer; x+=xVecBox){
				for(double y=-yVecBox*nLayer; y<=yVecBox*nLayer; y+=yVecBox){
					for(double z=-zVecBox*nLayer; z<=zVecBox*nLayer; z+=zVecBox){
						lsPosition.E(new double[]{x, y, z});
						translateBy.setTranslationVector(lsPosition);
						atomGroupActionTranslate.actionPerformed(molecule1);
	
						sum += potential.energy(pair);
						
						lsPosition.TE(-1);
						translateBy.setTranslationVector(lsPosition);
						atomGroupActionTranslate.actionPerformed(molecule1);
					}	
				}	
			}
		} else {
			sum += potential.energy(pair);
		}
		
		return sum;
		
	}
	
	public double d2phi_du2(int[] moleculei, int[] d) {
		
		/*
		 * d[i] = the i-th elements of the derivative
		 * u[i] = the generalized coordinate/ relative displacement from 
		 * 			nominal position 
		 */
        deltaU[0] = fixedDeltaU;
        deltaU[1] = fixedDeltaU;
        a[0][0] = computeA(moleculei, d, deltaU);
        
		double err = big;
        double d2fdu2 = Double.NaN;

		for(int i=1; i<ntab; i++){
			deltaU[0] = deltaU[0] /con;
			deltaU[1] = deltaU[1] /con;
			
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

	public double computeA(int[] moleculei, int[] d, double[] deltaU){
		
		if(d.length > 2){
			throw new RuntimeException("In computeA method; Derivative method can only do first and" +
					"second order derivative!!!");
		}
	
		/*
		 *  SECOND PARTIAL DERIVATIVE
		 *  
		 *  f is the energy function
		 */
		
		if(moleculei[0]==moleculei[1] && d[0]==d[1]){
			/*
			 * 
			 *  when u_i = u_j
			 *  
			 *      d2f            f(u_i + deltaU) - 2*f(u_i) + f(u_i - deltaU) 
			 *  -----------  = ----------------------------------------------------- 
			 *   du_i*du_i                            deltaU^2
			 *  
			 */
				
			double f_init = f(moleculei, generalizedCoord);
				
			generalizedCoord[0][d[0]] = deltaU[0];
			double f_up1 = f(moleculei, generalizedCoord);
	
			generalizedCoord[0][d[0]] = -deltaU[0];
			double f_um1 = f(moleculei, generalizedCoord);
				
			generalizedCoord[0][d[0]] = 0.0;
			generalizedCoord[1][d[1]] = 0.0;
				
			setToInitialPosition(moleculei);
			return (f_up1 - 2*f_init + f_um1)/(deltaU[0]*deltaU[1]);  
		
			
		} else {
				
			/*
			 *  when u_i, u_j
			 *  
			 *                  [  f(u_i+deltaU, u_j+deltaU) - f(u_i+deltaU, u_j-deltaU) - 
			 *     d2f                                    f(u_i-deltaU, u_j+deltaU) + f(u_i-deltaU, u_j-deltaU) ]
			 *  -----------  =  ------------------------------------------------------------------------------------ 
			 *   du_i*du_j                            4 * deltaU^2
			 *  
			 */
			
			generalizedCoord[0][d[0]] =  deltaU[0];
			generalizedCoord[1][d[1]] =  deltaU[1];
			double f_ip1jp1 = f(moleculei, generalizedCoord);

			generalizedCoord[0][d[0]] =  deltaU[0];
			generalizedCoord[1][d[1]] = -deltaU[1];
			double f_ip1jm1 = f(moleculei, generalizedCoord);

			generalizedCoord[0][d[0]] = -deltaU[0];
			generalizedCoord[1][d[1]] = -deltaU[1];
			double f_im1jm1 = f(moleculei, generalizedCoord);

			generalizedCoord[0][d[0]] = -deltaU[0];
			generalizedCoord[1][d[1]] =  deltaU[1];
			double f_im1jp1 = f(moleculei, generalizedCoord);

			generalizedCoord[0][d[0]] = 0.0;
			generalizedCoord[1][d[1]] = 0.0;
			
			setToInitialPosition(moleculei);
			return (f_ip1jp1 - f_ip1jm1 - f_im1jp1 + f_im1jm1)/(4*deltaU[0]*deltaU[1]);
		}
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

	protected CoordinateDefinitionNitrogen coordinateDefinition;
	protected P2Nitrogen potential;
	protected AtomActionTranslateBy translateBy;
	protected MoleculeChildAtomAction atomGroupActionTranslate;
	protected IVectorMutable lsPosition;
	protected double errt, fac, xVecBox, yVecBox, zVecBox;
	protected double[] deltaU = new double[2];
	protected double [][] a, generalizedCoord;
	protected boolean doLatticeSum = false;
	double fixedDeltaU = 0.01;
	final int ntab = 10;
	final double con = 1.4;
	final double con2 = con*con;
	final double big = Double.MAX_VALUE;
	final double safe = 2.0;
	
	
}
