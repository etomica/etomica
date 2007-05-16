package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisMonatomic;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveFcc;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2SoftSpherical;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.statmech.LennardJones;
import etomica.units.Energy;
import etomica.util.FunctionGeneral;

/**
 * Properties of a system of monatomic molecules occupying a lattice and interacting according
 * to a spherically-symmetric pair potential.  Properties are given by a lattice-dynamics treatment.
 * 
 * @author kofke
 *
 */
public class HarmonicCrystal {

    public HarmonicCrystal(int[] nCells, Primitive primitive, Basis basis, Potential2SoftSpherical potential) {
        this.potential = potential;
        this.nCells = (int[])nCells.clone();
        lattice = new BravaisLatticeCrystal(primitive, basis);
//        normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
        normalModes = new NormalModesPotential(nCells, primitive, basis, potential);
    }
    
    public double getLatticeEnergy() {
        FunctionGeneral function = new FunctionGeneral() {
            public Data f(Object obj) {
                data.x = potential.u(((Vector3D)obj).squared());
                return data;
            }
            public IDataInfo getDataInfo() {
                return dataInfo;
            }
            final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
            final DataDouble data = new DataDouble();
        };
        LatticeSumCrystal summer = new LatticeSumCrystal(lattice);
        summer.setMaxElement(49);
        summer.setK(lattice.getSpace().makeVector());
//        System.out.println("\n k:"+kVector.toString());
        double sum = 0;
        double basisDim = lattice.getBasis().getScaledCoordinates().length;
        DataGroupLSC data = (DataGroupLSC)summer.calculateSum(function);
        for(int j=0; j<basisDim; j++) {
            for(int jp=0; jp<basisDim; jp++) {
                sum += ((DataDouble)data.getDataReal(j,jp)).x; 
            }
        }
        return 0.5*sum/basisDim;
//            ((Tensor)sum[0]).map(chopper);
//            ((Tensor)sum[1]).map(chopper);
//            ((Tensor)sum[0]).ME(sum0);
//            System.out.println(sum[0].toString());
 //           System.out.println();
//            System.out.println(sum[1].toString());
    }
    
    public double getHelmholtzFreeEnergy(double temperature) {
        
        double coeffSum = 0.0;
//      int nA = (int)(2*coeffSum + 1);
        int D = lattice.getSpace().D();
        int cellCount = 1;
        int differ = 1;
        for(int i=0; i<D; i++) {
            cellCount *= nCells[i];
            if(nCells[i] % 2 == 0) {
                differ *= 2;
            }
        }
        int basisDim = lattice.getBasis().getScaledCoordinates().length;
        int moleculeCount = cellCount*basisDim;
        System.out.println("Space dimension: " + D);
        System.out.println("coeffSum: "+coeffSum);
        System.out.println("cellCount: "+cellCount);
        System.out.println("basisDim: "+basisDim);
        System.out.println("moleculeCount: "+moleculeCount);
        double jacobian = 0.5*D*(basisDim*(cellCount - differ)*Math.log(2.0) - Math.log(cellCount));
        System.out.println("differ, jacobian: " + differ + "\t" + jacobian);

        double[][] omega2 = normalModes.getOmegaSquared(null);//need to change signature of this method
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        
        double sumA = 0.0;
        double normalModeSum = 0.0;
        for(int k=0; k<omega2.length; k++) {
            double coeff = coeffs[k];
            for(int i=0; i<omega2[k].length; i++) {
                if(omega2[k][i] != 0.0) {
                    sumA += coeff*Math.log(omega2[k][i]*coeff/(temperature*Math.PI));
                    normalModeSum += coeff;
                }
            }
        }

        System.out.println("2*normalModeSum + D: " + (2*normalModeSum+D));
        System.out.println("D * moleculeCount: " + (D*moleculeCount));
        sumA -= jacobian;
        sumA /= moleculeCount;
        sumA *= temperature;
        sumA += getLatticeEnergy();
        return sumA;
    }
    
    public void setCellDensity(double newDensity) {
        double oldVolume = lattice.getPrimitive().unitCell().getVolume();
        double scale = newDensity * oldVolume;
        Primitive primitive = lattice.getPrimitive();
        primitive.scaleSize(1.0/Math.pow(scale, 1.0/lattice.getSpace().D()));
//        normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
        normalModes = new NormalModesPotential(nCells, primitive, lattice.getBasis(), potential);

    }
    
    public static void main(String[] args) {
        double T = 1;
        double rho = 1.0;
        PrimitiveFcc primitive = new PrimitiveFcc(Space3D.getInstance());
        Basis basis = new BasisMonatomic(Space3D.getInstance());

//        Primitive primitive = new PrimitiveCubic(Space3D.getInstance());
//        Basis basis = new BasisCubicFcc();
        
        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);

        int nC = 3;
        int[] nCells = new int[] {nC, nC, nC};
        
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential);
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
        
        System.out.println("Density: " + rho);
        System.out.println("Temperature: " + T);
        
        double u = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy(T);
        System.out.println("Lattice Energy: " + u);
        double uEos = LennardJones.uStaticFcc(rho);
        System.out.println("Energy from EOS: " + uEos);
        System.out.println("Helmholtz: " + a);
        double aEos = LennardJones.aResidualFcc(T,rho) + T*Math.log(rho) - 1.0;
        System.out.println("Helmholtz from EOS: " + aEos);
    }
    
    private NormalModes normalModes;
    private BravaisLatticeCrystal lattice;
    private int[] nCells;
    private Potential2SoftSpherical potential;
    private static final long serialVersionUID = 1L;
    
}
