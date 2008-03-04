package etomica.normalmode;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.IDataInfo;
import etomica.data.types.DataTensor;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.box.Box;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Tensor;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.units.Dimension;
import etomica.util.Arrays;
import etomica.util.Function;
import etomica.util.FunctionGeneral;
import etomica.util.RandomNumberGenerator;

/**
 * Uses analysis of 2nd derivatives to compute the normal modes for a Bravais lattice with a basis, 
 * occupied by atoms that interact with a simple, spherically-symmetric soft potential.
 */

public class NormalModesPotential implements NormalModes {

    public NormalModesPotential(int[] nCells, Primitive primitive, Basis basis, Potential2SoftSpherical potential, Space space) {
        
        harmonicFudge = 1.0;
        needToCalculateModes = true;
        
        lattice = new BravaisLatticeCrystal(primitive, basis);
        setPotential(potential);
        
        int nSites = nCells[0]*nCells[1]*nCells[2];
        Boundary boundary = new BoundaryDeformableLattice(primitive, new RandomNumberGenerator(), nCells);
        
        IBox box = new Box(boundary, space);

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
        //this function returns phi_{alpha,beta}, as defined in Dove Eq. 6.15
        FunctionGeneral function = new FunctionGeneral() {
            public Data f(Object obj) {
                Vector3D r = (Vector3D)obj;
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
            final DataTensor tensor = new DataTensor(Space3D.getInstance());
            final DataInfo dataInfo = new DataTensor.DataInfoTensor("",Dimension.MIXED,Space3D.getInstance());
            final Tensor3D identity = new Tensor3D(new double[] {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0});
        };
        
        LatticeSumCrystal summer = new LatticeSumCrystal(lattice);
        summer.setMaxLatticeShell(maxLatticeShell);
        IVector kVector = lattice.getSpace().makeVector();
        
        //calculation of self term
        kVector.E(0.0);
        summer.setK(kVector);
        System.out.println("\n k:"+kVector.toString()+"   in NormalModesPotential");
        DataGroupLSC sum0 = (DataGroupLSC)summer.calculateSum(function);
        Function chopper = new Function.Chop(1e-9);
        sum0.map(chopper);
//        System.out.println(sum0.toString());
        double[][] array = new double[eDim][eDim];
        for(int k=0; k<kDim; k++) {
            kVector.E(kFactory.getWaveVectors()[k]);
            summer.setK(kVector);
            System.out.println("k:"+kVector.toString());
            DataGroupLSC sum = (DataGroupLSC)summer.calculateSum(function);
            sum.map(chopper);
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                    sum.getDataReal(j,j).ME(sum0.getDataReal(j, jp));
                }
            }
            for(int j=0; j<basisDim; j++) {
                for(int jp=0; jp<basisDim; jp++) {
                    Tensor tensor = ((DataTensor)sum.getDataReal(j,jp)).x;
                    for(int alpha=0; alpha<spaceDim; alpha++) {
                        for(int beta=0; beta<spaceDim; beta++) {
                            array[spaceDim*j+alpha][spaceDim*jp+beta] = tensor.component(alpha, beta);
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
        
    public void setPotential(Potential2SoftSpherical potential) {
        needToCalculateModes = true;
        this.potential = potential;
    }

    public double[][][] getEigenvectors(IBox box) {
        if(needToCalculateModes) {
            calculateModes();
        }
        return eigenvectors;
    }

    public double[][] getOmegaSquared(IBox box) {
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

    private final BravaisLatticeCrystal lattice;
    private Potential2SoftSpherical potential;
    private WaveVectorFactory kFactory;
    private int maxLatticeShell;
    private double harmonicFudge;
    private double temperature;
    private double[][] omega2;
    private double[][][] eigenvectors;
    private boolean needToCalculateModes;

}
