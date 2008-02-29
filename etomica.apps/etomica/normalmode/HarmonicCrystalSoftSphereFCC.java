 package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.LatticeSumCrystal;
import etomica.lattice.LatticeSumCrystal.DataGroupLSC;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.potential.P2SoftSphere;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.Energy;
import etomica.util.Arrays;
import etomica.util.FunctionGeneral;

/**
 * Properties of a system of monatomic molecules occupying a lattice and interacting according
 * to a soft-sphere spherically-symmetric pair potential.  Properties are given by a lattice-dynamics treatment.
 * 
 * FCC Crystal Structure
 * 
 * @author kofke
 *
 */
public class HarmonicCrystalSoftSphereFCC {

    public HarmonicCrystalSoftSphereFCC(int[] nCells, Primitive primitive,
    		    Basis basis, Potential2SoftSpherical potential, Space _space) {
        this.potential = potential;
        this.nCells = (int[])nCells.clone();
        this.space = _space;
        lattice = new BravaisLatticeCrystal(primitive, basis);
        normalModes = new NormalModesPotential(nCells, primitive, basis, potential, space);
        setMaxLatticeShell(49);
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
        summer.setMaxLatticeShell(maxLatticeShell);
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
        System.out.println("cellCount: "+cellCount);
        System.out.println("basisDim: "+basisDim);
        System.out.println("moleculeCount: "+moleculeCount);
        double jacobian = 0.5*D*(basisDim*(cellCount - differ)*Math.log(2.0) - Math.log(cellCount));
        System.out.println("differ, jacobian: " + differ + "\t" + jacobian);

        double[][] omega2 = normalModes.getOmegaSquared(null);//need to change signature of this method
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        System.out.println("coeffs: "+Arrays.toString(coeffs));
        double sumA = 0.0;
        double normalModeSum = 0.0;
        double omega2zeroCount = 0;
        for(int k=0; k<omega2.length; k++) {
            double coeff = coeffs[k];
            for(int i=0; i<omega2[k].length; i++) {
                if(omega2[k][i] > 1.e-9) {
                    sumA += coeff*Math.log(omega2[k][i]*coeff/(temperature*Math.PI));
                    normalModeSum += coeff;
                } else {
                    omega2zeroCount++;
                }
            }
        }

        System.out.println("omega2==0 count: "+omega2zeroCount);
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
        normalModes = new NormalModesPotential(nCells, primitive, lattice.getBasis(), potential, space);
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }
    
    public int getMaxLatticeShell() {
        return maxLatticeShell;
    }

    public void setMaxLatticeShell(int maxLatticeShell) {
        this.maxLatticeShell = maxLatticeShell;
        normalModes.setMaxLatticeShell(maxLatticeShell);
    }

    public static void main(String[] args) {
        double rho = 1.338;
        double softness = 0.10;
        int maxLatticeShell = 49;
        int nC =3;
//        Primitive primitive = new PrimitiveFcc(Space3D.getInstance());
//        Basis basis = new BasisMonatomic(Space3D.getInstance());
        
        if (args.length > 0) {
            rho = Double.parseDouble(args[0]);
        }
        if (args.length > 1) {
            softness = Double.parseDouble(args[1]);
        }
        if (args.length > 2) {
            nC = Integer.parseInt(args[2]);
        }
        
        Primitive primitive = new PrimitiveCubic(Space3D.getInstance());
        Basis basis = new BasisCubicFcc();
        
        Space sp = Space3D.getInstance();
        final Potential2SoftSpherical potential = new P2SoftSphere(sp, 1.0, 1.0, softness);

        int[] nCells = new int[] {nC, nC, nC};
        
        HarmonicCrystalSoftSphereFCC harmonicCrystal = new HarmonicCrystalSoftSphereFCC(nCells, primitive, basis, potential, sp);
        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
        harmonicCrystal.setMaxLatticeShell(harmonicCrystal.maxLatticeShell);
        
        System.out.println("Density: " + rho);
        System.out.println("Softness: " + softness);
        
        double u = harmonicCrystal.getLatticeEnergy();
        
        double[] a = new double[100];
        double temp =0.01;
        for (int i =0; i<a.length; i++){
        	a[i] = harmonicCrystal.getHelmholtzFreeEnergy(temp);
        	temp+=0.01;
        }
        System.out.println(" ");
        
        temp=0.01;
        for (int i =0; i<a.length; i++){
        	System.out.println("Helmholtz Free Energy at T"+temp+ " is: "+a[i]);
        	temp+=0.01;
        }
        
        System.out.println("\nLattice Energy: " + u);

        
//        double latticeConstant = 1.0;
//        primitive = new PrimitiveHexagonal(Space3D.getInstance(), latticeConstant, Math.sqrt(8.0/3.0)*latticeConstant);
//        basis = new BasisHcp();
//        harmonicCrystal = new HarmonicCrystal(nCells, primitive, basis, potential);
//        harmonicCrystal.setCellDensity(rho/basis.getScaledCoordinates().length);
//        harmonicCrystal.setMaxLatticeShell(maxLatticeShell);
//        u = harmonicCrystal.getLatticeEnergy();
//        System.out.println("Lattice energy (HCP): "+u);
    }
    
    private NormalModesPotential normalModes;
    private BravaisLatticeCrystal lattice;
    private int[] nCells;
    private int maxLatticeShell;
    private Potential2SoftSpherical potential;
    private final Space space;
    private static final long serialVersionUID = 1L;
    
}
