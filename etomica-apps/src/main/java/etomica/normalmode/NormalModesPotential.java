/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataInfo;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.*;
import etomica.space3d.Tensor3D;
import etomica.units.dimensions.Dimension;
import etomica.math.function.Function;
import etomica.data.FunctionData;

/**
 * Uses analysis of 2nd derivatives to compute the normal modes for a Bravais lattice with a basis, 
 * occupied by atoms that interact with a simple, spherically-symmetric soft potential.
 */

public class NormalModesPotential implements NormalModes {

    public NormalModesPotential(int[] nCells, Primitive primitive, Basis basis, Potential2SoftSpherical potential, Space space) {
        this.space = space;
        needToCalculateModes = true;
        
        lattice = new BravaisLatticeCrystal(primitive, basis);
        setPotential(potential);
        
        int nSites = nCells[0]*nCells[1]*nCells[2];
        Boundary boundary = new BoundaryDeformableLattice(primitive, nCells);
        
        Box box = new Box(boundary, space);

//        System.out.println("Cell Density::: "+nSites/boundary.volume());
  //      System.out.println("Site Density: "+nSites/boundary.volume()*basis.getScaledCoordinates().length);
        kFactory = new WaveVectorFactorySimple(primitive, space);
        kFactory.makeWaveVectors(box);
  //      System.out.println("Number of wave vectors: "+kFactory.getWaveVectors().length);
        
        double sum = 0.0;
        for(int i=0; i<kFactory.getWaveVectors().length; i++) {
            sum += kFactory.getCoefficients()[i];
        }
   //     System.out.println("nCells: "+ Arrays.toString(nCells));
     //   System.out.println("Number of wave vectors represented: "+2.0*sum);
    }
    
    public void calculateModes() {
        int kDim = kFactory.getWaveVectors().length;
        int spaceDim = lattice.getSpace().D();
        int basisDim = lattice.getBasis().getScaledCoordinates().length;
        int eDim = basisDim * spaceDim;
        omega2 = new double[kDim][eDim];
        eigenvectors = new double[kDim][eDim][eDim];
        //this function returns phi_{alpha,beta}, as defined in Dove Eq. 6.15
        FunctionData<Object> function = new FunctionData<Object>() {
            public IData f(Object obj) {
                Vector r = (Vector)obj;
                tensor.x.Ev1v2(r, r);
                double r2 = r.squared();
                double dW = potential.du(r2);
                double d2W = potential.d2u(r2);
                tensor.TE(1.0/(r2*r2)*(dW - d2W));
                tensor.x.PEa1Tt1(-dW/r2,identity);
                return tensor;
            }
            public IDataInfo getDataInfo() {
                return dataInfo;
            }
            final DataTensor tensor = new DataTensor(space);
            final DataInfo dataInfo = new DataTensor.DataInfoTensor("", Dimension.MIXED, space);
            final Tensor identity = new Tensor3D(new double[][] {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}});
        };
        
        LatticeSumCrystal summer = new LatticeSumCrystal(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        Vector kVector = lattice.getSpace().makeVector();
        
        
        
   
        	
        //calculation of self term
        kVector.E(0.0);
        summer.setK(kVector);
 //       System.out.println("\n k:"+kVector.toString()+"   in NormalModesPotential");
        DataGroupLSC sum0 = (DataGroupLSC)summer.calculateSum(function);
        Function chopper = new Function.Chop(1e-9);
        sum0.map(chopper);
//        System.out.println(sum0.toString());
        double[][] array = new double[eDim][eDim];
        
        try{

        	FileWriter fileWriterK = null;
        	FileWriter fileWriterVal = null;
        	FileWriter fileWriterVec = null;

        	if (fileName != null) {
        	    fileWriterK = new FileWriter(getFileName()+".k");
                fileWriterVal = new FileWriter(getFileName()+".val");
                fileWriterVec = new FileWriter(getFileName()+".vec");
        	}
        	
        	
        double[] kCoefficients = kFactory.getCoefficients();	
        	
        for(int k=0; k<kDim; k++) {
            kVector.E(kFactory.getWaveVectors()[k]);
            
            // output .k file
            if (fileName != null) {
                fileWriterK.write(Double.toString(kCoefficients[k]));
                for (int n=0; n< kFactory.getWaveVectors()[k].getD(); n++){
                	fileWriterK.write(" "+ kFactory.getWaveVectors()[k].getX(n));
                }
                fileWriterK.write("\n");
            }
            
            
            summer.setK(kVector);
  //          System.out.println("k:"+kVector.toString());
            DataGroupLSC sum = (DataGroupLSC)summer.calculateSum(function);
            sum.map(chopper);
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                    sum.getDataReal(j,j).ME(sum0.getDataReal(j, jp));
                }
            }
            // impose symmetry on the matrix.  elements of the tensor can have elements that are different
            // only in the last digit (due numerical precision issues).  With a not exactly symmetric
            // matrix Jama assumes it's asymmetric and finds not-orthogonal eigenvectors.
            for(int j=0; j<basisDim; j++) {
                for(int jp=j+1; jp<basisDim; jp++) {
                    // grab mirror blocks for j,jp and jp,j
                    Tensor tensorj_jp = ((DataTensor)sum.getDataReal(j,jp)).x;
                    Tensor tensorjp_j = ((DataTensor)sum.getDataReal(j,jp)).x;
                    for(int alpha=0; alpha<spaceDim; alpha++) {
                        for(int beta=0; beta<spaceDim; beta++) {
                            // average opposite components from opposite blocks
                            double v = 0.5*(tensorj_jp.component(alpha, beta) + tensorjp_j.component(beta, alpha));
                            array[spaceDim*j+alpha][spaceDim*jp+beta] = v;
                            array[spaceDim*jp+beta][spaceDim*j+alpha] = v;
                        }
                    }
                }
                // grab diagonal block
                Tensor tensor = ((DataTensor)sum.getDataReal(j,j)).x;
                for(int alpha=0; alpha<spaceDim; alpha++) {
                    for(int beta=alpha+1; beta<spaceDim; beta++) {
                        // average opposite components
                        double v = 0.5*(tensor.component(alpha, beta) + tensor.component(beta, alpha));
                        array[spaceDim*j+alpha][spaceDim*j+beta] = v;
                        array[spaceDim*j+beta][spaceDim*j+alpha] = v;
                    }
                    // grab diagonal component
                    array[spaceDim*j+alpha][spaceDim*j+alpha] = tensor.component(alpha, alpha);
                }
            }

            Matrix matrix = new Matrix(array);
            EigenvalueDecomposition ed = matrix.eig();
            double[] eVals = ed.getRealEigenvalues();
            double[][] eVecs = ed.getV().getArray();
            
  //          System.out.println("Real eigenvalues: " + Arrays.toString(eVals));
            
            if (fileName != null) {
                // output .val file
                for (int ival=0; ival<eVals.length; ival++){
                	if (eVals[ival] < 1E-10){
                		fileWriterVal.write("0.0 ");
                	} else {
                		fileWriterVal.write(1/eVals[ival]+ " ");
                	}
                }
                fileWriterVal.write("\n");

                // output .vec file
                for (int ivec=0; ivec<eDim; ivec++ ){
                	for(int jvec=0; jvec<eDim; jvec++){
                		if (Math.abs(eVecs[jvec][ivec])<1e-10){
                			fileWriterVec.write("0.0 ");
                		} else {
                			fileWriterVec.write(eVecs[jvec][ivec] + " ");
                		}
                	}
                	fileWriterVec.write("\n");
                }
            }
            
//            System.out.println("Imag eigenvalues: " + Arrays.toString(ed.getImagEigenvalues()));
            
            for(int j=0; j<eDim; j++) {
                if (eVals[j] < 1E-12) {
                    omega2[k][j] = Double.POSITIVE_INFINITY;
                }
                else {
                    omega2[k][j] = eVals[j];
                }
                for(int m=0; m<eDim; m++) {
                    eigenvectors[k][j][m] = eVecs[m][j];//need to check if indexes are right
                }
            }
//            System.out.println(sum[0].toString());
//            System.out.println();
//            System.out.println(sum[1].toString());
        }
        if (fileName != null) {
            fileWriterK.close();
            fileWriterVal.close();
            fileWriterVec.close();
        }
    }
        catch (IOException e) {
            throw new RuntimeException(e);
        }
        
        needToCalculateModes = false;
    }
        
    public void setPotential(Potential2SoftSpherical potential) {
        needToCalculateModes = true;
        this.potential = potential;
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
        // we ignore fudge
    }
    
    public void setTemperature(double newTemperature) {
        // we ignore temperature
    }
    
	public String getFileName() {
		return fileName;
	}

	public void setFileName(String filename) {
		this.fileName = filename;
	}

    protected final Space space;
    private final BravaisLatticeCrystal lattice;
    private Potential2SoftSpherical potential;
    private WaveVectorFactory kFactory;
    private int maxLatticeShell;
    private double[][] omega2;
    private double[][][] eigenvectors;
    private boolean needToCalculateModes;
    private String fileName;
}
