package etomica.normalmode;

import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.IDataInfo;
import etomica.data.types.DataDouble;
import etomica.lattice.BravaisLattice;
import etomica.lattice.LatticeSum;
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

    public HarmonicCrystal(int[] nCells, Primitive primitive, Potential2SoftSpherical potential) {
        this.potential = potential;
        this.nCells = (int[])nCells.clone();
        lattice = new BravaisLattice(primitive);
        normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
    }
    
    public double getLatticeEnergy() {
        FunctionGeneral function = new FunctionGeneral() {
            public Data f(Object obj) {
                Vector3D r = (Vector3D)obj;
                double r2 = r.squared();
                data.x = potential.u(r2);
                return data;
            }
            public IDataInfo getDataInfo() {
                return dataInfo;
            }
            final DataInfo dataInfo = new DataDouble.DataInfoDouble("Lattice energy", Energy.DIMENSION);
            final DataDouble data = new DataDouble();
        };
        LatticeSum summer = new LatticeSum(lattice);
        summer.setMaxElement(49);
        summer.setK(lattice.getSpace().makeVector());
//        System.out.println("\n k:"+kVector.toString());
        return 0.5*((DataDouble)summer.calculateSum(function).getData(0)).x;
//            ((Tensor)sum[0]).map(chopper);
//            ((Tensor)sum[1]).map(chopper);
//            ((Tensor)sum[0]).ME(sum0);
//            System.out.println(sum[0].toString());
 //           System.out.println();
//            System.out.println(sum[1].toString());
    }
    
    public double getHelmholtzFreeEnergy(double temperature) {
        
        double sumA = 0.0;
        
        double[][] omega2 = normalModes.getOmegaSquared(null);//need to change signature of this method
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        double coeffSum = 0.0;
        for(int k=0; k<omega2.length; k++) {
            double coeff = coeffs[k];
            coeffSum += coeff;
            for(int i=0; i<omega2[k].length; i++) {
                if(omega2[k][i] != 0.0) {
                    sumA += coeff*Math.log(omega2[k][i]*coeff/(temperature*Math.PI));
                }
            }
        }
        int nA = (int)(2*coeffSum + 1);
        System.out.println("nA: "+nA);
        int D = lattice.getSpace().D();
        double AHarmonic = sumA;
        if (nA % 2 == 0) {
            if (D == 3) {
                AHarmonic += Math.log(((3*nA + 3)/2.0) / Math.pow(nA,1.5));
                //0.5*D*Math.log(nA) - 0.5*(nA-1)*Math.log(2.0*Math.PI)
            }
        }
        else {
            if (D == 3) {
                AHarmonic += Math.log(((3*nA - 18)/2.0) / Math.pow(nA,1.5));
            }
        }
//        if(nA % 2 == 0) AHarmonic += 0.5*Math.log(2.0);
        AHarmonic /= nA;
        AHarmonic *= temperature;
        AHarmonic += getLatticeEnergy();
        return AHarmonic;
    }
    
    public void setCellDensity(double newDensity) {
        double oldVolume = lattice.getPrimitive().unitCell().getVolume();
        double scale = newDensity * oldVolume;
        Primitive primitive = lattice.getPrimitive();
        primitive.scaleSize(1.0/Math.pow(scale, 1.0/lattice.getSpace().D()));
        normalModes = new NormalModesSoftSpherical(nCells, primitive, potential);
    }
    
    public static void main(String[] args) {
        double T = 1;
        double rho = 1.0;
        PrimitiveFcc primitive = new PrimitiveFcc(Space3D.getInstance());
        primitive.setCubicSize(Math.pow(2.0, 1.0/6.0));//unit density
        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);

        int nC = 4;
        int[] nCells = new int[] {nC, nC, nC};
        
        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(nCells, primitive, potential);
        harmonicCrystal.setCellDensity(rho);
        
        System.out.println("Density: " + rho);
        System.out.println("Temperature: " + T);
        
        double u = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy(T);
        System.out.println("Lattice Energy: " + u);
        System.out.println("Helmholtz: " + a);
        double aEos = LennardJones.aResidualFcc(T,rho) + T*Math.log(rho) - 1.0;
        double uEos = LennardJones.uStaticFcc(rho);
        System.out.println("Energy from EOS: " + uEos);
        System.out.println("Helmholtz from EOS: " + aEos);
    }
    
    private NormalModes normalModes;
    private BravaisLattice lattice;
    private int[] nCells;
    private Potential2SoftSpherical potential;
    private static final long serialVersionUID = 1L;
    
}
