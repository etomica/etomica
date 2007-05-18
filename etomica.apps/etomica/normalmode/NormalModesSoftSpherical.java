package etomica.normalmode;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.IDataInfo;
import etomica.data.types.DataGroup;
import etomica.data.types.DataTensor;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeSum;
import etomica.lattice.crystal.Primitive;
import etomica.phase.Phase;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Boundary;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.IVector;
import etomica.space3d.Space3D;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.units.Dimension;
import etomica.util.Arrays;
import etomica.util.Function;
import etomica.util.FunctionGeneral;
import etomica.util.RandomNumberGenerator;

/**
 * Computes the normal modes for a Bravais lattice occupied by atoms that interact
 * with a simple, spherically-symmetric soft potential, using analysis of second
 * derivatives.
 */

public class NormalModesSoftSpherical implements NormalModes {

    public NormalModesSoftSpherical(int[] nCells, Primitive primitive, Potential2SoftSpherical potential) {
        
        harmonicFudge = 1.0;
        needToCalculateModes = true;
        
        setPrimitive(primitive);
        setPotential(potential);
        
        int nSites = nCells[0]*nCells[1]*nCells[2];
        Boundary boundary = new BoundaryDeformableLattice(primitive, new RandomNumberGenerator(), nCells);
        
        Phase phase = new Phase(boundary);

        System.out.println("Density: "+nSites/boundary.volume());
        kFactory = new WaveVectorFactorySimple(primitive);
        kFactory.makeWaveVectors(phase);
        System.out.println("Number of wave vectors: "+kFactory.getWaveVectors().length);
        
        double sum = 0.0;
        for(int i=0; i<kFactory.getWaveVectors().length; i++) {
            sum += kFactory.getCoefficients()[i];
        }
        System.out.println("nCells: "+Arrays.toString(nCells));
        System.out.println("Number of wave vectors represented: "+2.0*sum);
    }
    
    public void calculateModes() {
        omega2 = new double[kFactory.getWaveVectors().length][lattice.D()];
        eigenvectors = new double[omega2.length][lattice.D()][lattice.D()];
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
        
        LatticeSum summer = new LatticeSum(lattice);
        summer.setMaxElement(49);
        IVector kVector = new Vector3D();
        kVector.E(0.0);
        summer.setK(kVector);
        System.out.println("\n k:"+kVector.toString()+"   in NormalModeSoftSpherical");
        DataGroup sum = summer.calculateSum(function);
        DataTensor sum0 = new DataTensor(Space3D.getInstance());
        sum0.E(sum.getData(0));
        Function chopper = new Function.Chop(1e-9);
        sum0.map(chopper);
//        System.out.println(sum0.toString());
        double[][] array = new double[3][3];
        for(int k=kFactory.getWaveVectors().length-1; k>0; k--) {
            kVector.E(kFactory.getWaveVectors()[k]);
            summer.setK(kVector);
            System.out.println("k:"+kVector.toString());
            sum = summer.calculateSum(function);
            sum.map(chopper);
            sum.getData(0).ME(sum0);
            ((DataTensor)sum.getData(0)).x.assignTo(array);
            Matrix matrix = new Matrix(array);
            EigenvalueDecomposition ed = matrix.eig();
            double[] eVals = ed.getRealEigenvalues();
            double[][] eVecs = ed.getV().getArray();
            System.out.println("Real eigenvalues: " + Arrays.toString(eVals));
//            System.out.println("Imag eigenvalues: " + Arrays.toString(ed.getImagEigenvalues()));
            
            for(int j=0; j<eVals.length; j++) {
                omega2[k][j] = eVals[j];
                for(int m=0; m<eVals.length; m++) {
                    eigenvectors[k][j][m] = eVecs[j][m];//need to check if indexes are right
                }
            }
//            System.out.println(sum[0].toString());
//            System.out.println();
//            System.out.println(sum[1].toString());
        }
        needToCalculateModes = false;
    }
    
    public void setPrimitive(Primitive primitive) {
        needToCalculateModes = true;
        lattice = new BravaisLattice(primitive);
    }
    
    public void setPotential(Potential2SoftSpherical potential) {
        needToCalculateModes = true;
        this.potential = potential;
    }

    public double[][][] getEigenvectors(Phase phase) {
        if(needToCalculateModes) {
            calculateModes();
        }
        return eigenvectors;
    }

    public double[][] getOmegaSquared(Phase phase) {
        if(needToCalculateModes) {
            calculateModes();
        }
        return omega2;
    }

    public WaveVectorFactory getWaveVectorFactory() {
        return kFactory;
    }

    public void setHarmonicFudge(double newHarmonicFudge) {
        needToCalculateModes = true;
        harmonicFudge = newHarmonicFudge;
    }
    
    public void setTemperature(double newTemperature) {
        needToCalculateModes = true;
        temperature = newTemperature;
    }

    private BravaisLattice lattice;
    private Potential2SoftSpherical potential;
    private WaveVectorFactory kFactory;
    private double harmonicFudge;
    private double temperature;
    private double[][] omega2;
    private double[][][] eigenvectors;
    private boolean needToCalculateModes;

}
