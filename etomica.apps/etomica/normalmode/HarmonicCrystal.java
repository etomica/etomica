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

public class HarmonicCrystal {

    public HarmonicCrystal(Primitive primitive, Potential2SoftSpherical potential) {
        this.potential = potential;
        lattice = new BravaisLattice(primitive);
        normalModes = new NormalModesSoftSpherical(primitive, potential);
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
        summer.setK(lattice.getSpace().makeVector());
//        System.out.println("\n k:"+kVector.toString());
        return 0.5*((DataDouble)summer.calculateSum(function)[0]).x;
//            ((Tensor)sum[0]).map(chopper);
//            ((Tensor)sum[1]).map(chopper);
//            ((Tensor)sum[0]).ME(sum0);
//            System.out.println(sum[0].toString());
 //           System.out.println();
//            System.out.println(sum[1].toString());
    }
    
    public double getHelmholtzFreeEnergy() {
        
        double sumA = 0.0;
        
        double[][] omega2 = normalModes.getOmegaSquared(null);//need to change signature of this method
        double[] coeffs = normalModes.getWaveVectorFactory().getCoefficients();
        double coeffSum = 0.0;
        for(int k=0; k<omega2.length; k++) {
            double coeff = coeffs[k];
            coeffSum += coeff;
            for(int i=0; i<omega2[k].length; i++) {
                sumA += coeff*Math.log(omega2[k][i]*coeff);
            }
        }
        int nA = (int)(2*coeffSum + 1);
        double AHarmonic = 0.5*Math.log(nA) - 0.5*(nA-1)*Math.log(2.0*Math.PI) + sumA;
//        if(nA % 2 == 0) AHarmonic += 0.5*Math.log(2.0);
        AHarmonic /= nA;
        AHarmonic += getLatticeEnergy();
        System.out.println("Free energy: "+AHarmonic);
        return AHarmonic;
    }
    
    public static void main(String[] args) {
        double T = 1.0;
        double rho = 1.0;
        PrimitiveFcc primitive = new PrimitiveFcc(Space3D.getInstance());
        primitive.setCubicSize(Math.pow(2.0, 1.0/6.0));//unit density
        final Potential2SoftSpherical potential = new P2LennardJones(Space3D.getInstance(), 1.0, 1.0);

        HarmonicCrystal harmonicCrystal = new HarmonicCrystal(primitive, potential);
        
        double u = harmonicCrystal.getLatticeEnergy();
        double a = harmonicCrystal.getHelmholtzFreeEnergy();
        System.out.println("Lattice Energy: " + u);
        System.out.println("Helmholtz: " + a);
        double aEos = LennardJones.aFcc(T,rho);
        System.out.println("Helmholtz from EOS: " + aEos);
    }
    
    private NormalModes normalModes;
    private BravaisLattice lattice;
    private Potential2SoftSpherical potential;
    private static final long serialVersionUID = 1L;
    
}
