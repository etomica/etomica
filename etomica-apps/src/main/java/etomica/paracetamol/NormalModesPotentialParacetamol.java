/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.normalmode.NormalModes;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.paracetamol.LatticeSumCrystalParacetamol.DataGroupLSCParacetamol;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.units.dimensions.Dimension;
import etomica.util.Arrays;
import etomica.math.function.Function;
import etomica.math.numerical.FiniteDifferenceDerivative;

/**
 * Uses analysis of 2nd derivatives to compute the normal modes for a Bravais lattice with a basis, 
 * occupied by atoms that interact with a simple, spherically-symmetric soft potential.
 *
 * @author Tai Tan
 */

public class NormalModesPotentialParacetamol implements NormalModes {

    public NormalModesPotentialParacetamol(int[] nCells, Primitive primitive,
    		               Basis basis, Space space) {
        
        harmonicFudge = 1.0;
        needToCalculateModes = true;
        
        lattice = new BravaisLatticeCrystal(primitive, basis);
        
        int nSites = nCells[0]*nCells[1]*nCells[2];
        Boundary boundary = new BoundaryDeformableLattice(primitive, nCells);
        
        box = new Box(boundary, space);

        System.out.println("Cell Density: "+nSites/boundary.volume());
        System.out.println("Site Density: "+nSites/boundary.volume()*basis.getScaledCoordinates().length);
        kFactory = new WaveVectorFactorySimple(primitive, space);
        kFactory.makeWaveVectors(box);
        System.out.println("Number of wave vectors: "+kFactory.getWaveVectors().length);
        
        double sum = 0.0;
        for(int i=0; i<kFactory.getWaveVectors().length; i++) {
            sum += kFactory.getCoefficients()[i];
        }
        System.out.println("nCells: "+Arrays.toString(nCells));
        System.out.println("Number of wave vectors represented: "+2.0*sum);
    }
    
    public void calculateModes() {
    	
        int kDim = kFactory.getWaveVectors().length;
        int spaceDim = lattice.getSpace().D();
        int basisDim = lattice.getBasis().getScaledCoordinates().length;
        int eDim = basisDim * spaceDim;
        omega2 = new double[kDim][eDim];
        eigenvectors = new double[kDim][eDim][eDim];
        NumericalDerivativeEnergyParacetamol dW = new NumericalDerivativeEnergyParacetamol(box, potentialMaster);
       
        
        /*
         * this function returns phi_{alpha,beta}, as defined in Dove Eq. 6.15
         * alpha and beta each has 6 components
         * 
         * lattice2ndDerivative is the inner class that computes the 2nd derivatives
         * 	for the 6x6 matrix
         */
        Lattice2ndDerivative lattice2ndDerivative = new Lattice2ndDerivative(coordinateDefinitionParacetamol, dW);
        
        LatticeSumCrystalParacetamol summer = new LatticeSumCrystalParacetamol(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        Vector kVector = lattice.getSpace().makeVector();
        
        //calculation of self term
        kVector.E(0.0);
        summer.setK(kVector);
        System.out.println("\n k:"+kVector.toString()+"   in NormalModesPotentialParacetamol");
        summer.setCoordinateDefinitionParacetamol(coordinateDefinitionParacetamol);
        
        DataGroupLSCParacetamol sum0 = (DataGroupLSCParacetamol)summer.calculateSum(lattice2ndDerivative);
        
        System.out.println("NormalModesPotentialParacetamol!!");
        Function chopper = new Function.Chop(1e-9);
        sum0.map(chopper);
//      a System.out.println(sum0.toString());
        
        System.out.println(" after!");
        
        double[][] array = new double[eDim][eDim];
        for(int k=0; k<kDim; k++) {
            kVector.E(kFactory.getWaveVectors()[k]);
            summer.setK(kVector);
            System.out.println("k:"+kVector.toString());
            DataGroupLSCParacetamol sum = (DataGroupLSCParacetamol)summer.calculateSum(lattice2ndDerivative);
            sum.map(chopper);
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                    sum.getDataReal(j,j).ME(sum0.getDataReal(j, jp));
                }
            }
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                	
                	
                	/////////////////////////
                	DataDoubleArray tensor = ((DataDoubleArray)sum.getDataReal(j, jp));
                    //Tensor tensor = ((DataTensor)sum.getDataReal(j,jp)).x;
                    
                    
                    for(int alpha=0; alpha<spaceDim; alpha++) {
                        for(int beta=0; beta<spaceDim; beta++) {
                        	
                            array[spaceDim*j+alpha][spaceDim*jp+beta] = tensor.getValue(new int[]{alpha, beta});
                        }
                    }
                }
            }
            Matrix matrix = new Matrix(array);
            EigenvalueDecomposition ed = matrix.eig();
            double[] eVals = ed.getRealEigenvalues();
            double[][] eVecs = ed.getV().getArray();
            System.out.println("Real eigenvalues: " + Arrays.toString(eVals));
//            System.out.println("Imag eigenvalues: " + Arrays.toString(ed.getImagEigenvalues()));
            
            for(int j=0; j<eDim; j++) {
                omega2[k][j] = eVals[j];
                for(int m=0; m<eDim; m++) {
                    eigenvectors[k][j][m] = eVecs[m][j];//need to check if indexes are right
                }
            }
//            System.out.println(sum[0].toString());
//            System.out.println();
//            System.out.println(sum[1].toString());
        }
        needToCalculateModes = false;
    }
 


    public double[][][] getEigenvectors() {
        if(needToCalculateModes) {
            calculateModes();
        }
        return eigenvectors;
    }

    public double[][] getOmegaSquared() {
        if(needToCalculateModes) {
            calculateModes();
        }
        return omega2;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return kFactory;
    }

    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    public void setMaxLatticeShell(int maxLatticeShell) {
        this.maxLatticeShell = maxLatticeShell;
    }

    public void setHarmonicFudge(double newHarmonicFudge) {
        needToCalculateModes = true;
        harmonicFudge = newHarmonicFudge;
    }
    
    public void setTemperature(double newTemperature) {
        needToCalculateModes = true;
        temperature = newTemperature;
    }

    public void setCoordinateDefinition(CoordinateDefinitionParacetamol coordinateDefinitionParacetamol ){
    	this.coordinateDefinitionParacetamol = coordinateDefinitionParacetamol;
    }
    
    public CoordinateDefinitionParacetamol getCoordinateDefinition(){
    	return coordinateDefinitionParacetamol;
    }
    
	public void setPotentialMaster(PotentialMaster potentialMaster) {
		this.potentialMaster = potentialMaster;
	}
    
	public Box getBox() {
		return box;
	}

	public void setBox(Box box) {
		this.box = box;
	}
    
	public PotentialMaster getPotentialMaster() {
		return potentialMaster;
	}
	
    public static final class Lattice2ndDerivative implements LatticeEnergyParacetamol {
    	
    	public Lattice2ndDerivative(CoordinateDefinitionParacetamol coordinateDefinitionParacetamol, NumericalDerivativeEnergyParacetamol dW){
    		this.coordinateDefinitionParacetamol = coordinateDefinitionParacetamol;
    		this.dW = dW;
    		singleMolDim = coordinateDefinitionParacetamol.getCoordinateDim()
  		  	/coordinateDefinitionParacetamol.getBasisCells()[0].molecules.getMoleculeCount();
    		tag = new DataTag();
    		tensorDoubleArray = new DataDoubleArray(new int[]{singleMolDim, singleMolDim});
    		dataInfo = new DataDoubleArray.DataInfoDoubleArray("",Dimension.MIXED, new int[]{singleMolDim, singleMolDim});
  	
    	}
    	

		public IData getData() {
			//Object is passed in as a set of generalized coordinates within a cell
			double[] u = new double[coordinateDefinitionParacetamol.getCoordinateDim()];
			
			int[] d2 = new int[coordinateDefinitionParacetamol.getCoordinateDim()];
			dW.setCoordinateDefinition(coordinateDefinitionParacetamol);
		
			for (int i=indexjp*singleMolDim; i< (indexjp+1)*singleMolDim; i++){
				
				/*
				 * assignment of the first derivative
				 */
				d2[i]++;
				double[] columnArrayi = new double[singleMolDim];
				
				for (int k=indexj*singleMolDim; k<(indexj+1)*singleMolDim; k++){
					/*
		    		 * assignment of the second derivative
		    		 */
					d2[k]++;
					
		        	FiniteDifferenceDerivative d2W = new FiniteDifferenceDerivative(dW);
		    		d2W.setH(0.00001);
		    		
		    		//System.out.println("dimension of u is: "+ u.getLength());
		    		columnArrayi[k] = d2W.df(d2, u);
		    		d2[k]--;
		    	}
				
				System.out.println("NormalModesPotentialParacetamol secondDerivatives!!");
				
				tensorDoubleArray.assignColumnFrom(i, columnArrayi);
				d2[i]--;
			}
			return tensorDoubleArray;
         }
		
		/* (non-Javadoc)
		 * @see etomica.paracetamol.LatticeEnergyParacetamol#setMolecule(int, int)
		 */
		public void setMolecule(int indexj, int indexjp){
			this.indexj = indexj;
			this.indexjp = indexjp;
		}

		public IEtomicaDataInfo getDataInfo() {
			 return dataInfo;
		}

		public DataTag getTag(){
			return tag;
		}
		
		private final CoordinateDefinitionParacetamol coordinateDefinitionParacetamol;
		private final DataTag tag;
		private final DataDoubleArray tensorDoubleArray;
		private final DataInfo dataInfo;
		private final int singleMolDim;
		private final NumericalDerivativeEnergyParacetamol dW;
		private int indexj, indexjp;
	}

	
    private CoordinateDefinitionParacetamol coordinateDefinitionParacetamol;
    private PotentialMaster potentialMaster;
    private Box box;
    private final BravaisLatticeCrystal lattice;
    private WaveVectorFactory kFactory;
    private int maxLatticeShell;
    private double harmonicFudge;
    private double temperature;
    private double[][] omega2;
    private double[][][] eigenvectors;
    private boolean needToCalculateModes;

}
