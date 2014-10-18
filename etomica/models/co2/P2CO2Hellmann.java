package etomica.models.co2;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Oxygen;
import etomica.potential.IPotentialTorque;
import etomica.simulation.Simulation;
import etomica.space.ISpace;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.CompoundUnit;
import etomica.units.Kelvin;
import etomica.units.Null;
import etomica.units.Unit;

/**
 * Ab initio potential for CO2-CO2 developed by R. Hellmann.
 *
 * http://dx.doi.org/10.1016/j.cplett.2014.08.057
 * 
 * @author Andrew Schultz
 */
public class P2CO2Hellmann implements IPotentialTorque {

    protected static final double[] posA = new double[]{-1.28815171291, -1.17769797231, -0.18133162098, 0.00000000000, 0.18133162098, 1.17769797231, 1.28815171291};
    protected static final double[] posB = new double[]{-1.28741781626, -1.18192825424, -0.18607849166, 0.00000000000, 0.18607849166, 1.18192825424, 1.28741781626};
    protected static final double[] qA = new double[]{-191.214602050, 163.370217205, -2747.80131789, 5551.29140547, -2747.80131789, 163.370217205, -191.214602050};
    protected static final double[] qB = new double[]{-197.417207828, 168.070083318, -2559.64083227, 5177.97591356, -2559.64083227, 168.070083318, -197.417207828};
    protected static final int[] siteID = new int[]{0,1,2,3,2,1,0};
    protected static final double[] AA = new double[]{-0.206097356066E+07, 0.636857652220E+07, -0.194431965659E+08, 0.379428245681E+08, -0.112075555946E+07, 0.460535048972E+08, -0.118847782748E+09, -0.130239550952E+08, 0.425605973419E+08, -0.114999594389E+09};
    protected static final double[] AB = new double[]{-0.247910365353E+07, 0.659160470472E+07, -0.197776308389E+08, 0.384165630648E+08, -0.124570324466E+07, 0.451317323034E+08, -0.116048612008E+09, -0.103079402689E+08, 0.340824968085E+08, -0.915027698701E+08};
    protected static final double[] alphaA = new double[]{2.15463968686, 3.20330234880, 2.44546409876, 2.47207342773, 1.67720522533, 2.65642378772, 2.76521425413, 2.99122464510, 2.82606075114, 2.94148736356};
    protected static final double[] alphaB = new double[]{2.08319218048, 3.16681447768, 2.46163539534, 2.48589087370, 1.67813668662, 2.65969570294, 2.77169644514, 2.98535796569, 2.75870881239, 2.87267355769};
    protected static final double[] bA = new double[]{3.10392728450, 2.52351588006, 1.60598918542, 1.91489841285, 2.06812316859, 1.49177611494, 4.13165202187, 2.77767103339, 2.48867537323, 2.30146553200};
    protected static final double[] bB = new double[]{3.14980106637, 2.46903752251, 1.57103563097, 1.89845841233, 2.14451960163, 1.46843191121, 4.14021127755, 2.72634741238, 2.44815795987, 2.27614875317};
    protected static final double[] C6A = new double[]{-0.259827667988E+08, 0.606040735312E+08, -0.128063652046E+09, 0.106278254027E+09, -0.957062410858E+08, 0.729091325685E+08, 0.259781327401E+08, 0.146735186656E+09, -0.345229777876E+09, 0.709212843263E+09};
    protected static final double[] C6B = new double[]{-0.306747626563E+08, 0.698469835305E+08, -0.143806191593E+09, 0.121226824365E+09, -0.109398472925E+09, 0.811702881095E+08, 0.263241896284E+08, 0.126349448908E+09, -0.285769208067E+09, 0.551179708953E+09};
    protected static final double[] C8A = new double[]{0.189027654711E+09, -0.687384314073E+09, 0.312617758706E+10, -0.253453395293E+10, 0.930718903055E+09, -0.183474802396E+10, -0.175211033214E+09, -0.698068389813E+09, 0.191961111413E+10, -0.337571242220E+10};
    protected static final double[] C8B = new double[]{0.211522217149E+09, -0.810638994730E+09, 0.355929066714E+10, -0.286891373977E+10, 0.114677667224E+10, -0.210805303525E+10, -0.173569859005E+09, -0.496759975158E+09, 0.122323855871E+10, -0.131218053988E+10};
    protected static final int nsites = 7;
    
    protected final IVectorMutable ri, rj;
    protected final double[] pos, q;
    protected final double[][] A, alpha, b, C6, C8;
    
    public enum Parameters {
        A, B
    }
    
    public P2CO2Hellmann(ISpace space, Parameters param) {
        A = new double[4][4];
        alpha = new double[4][4];
        b = new double[4][4];
        C6 = new double[4][4];
        C8 = new double[4][4];
        q = new double[4];
        Unit sqrtK = new CompoundUnit(new Unit[]{Kelvin.UNIT}, new double[]{0.5});
        if (param == Parameters.A) {
            pos = posA;
            for (int i=0; i<q.length; i++) {
                q[i] = sqrtK.toSim(qA[i]);
            }
            ijInit(AA, A, Kelvin.UNIT);
            ijInit(alphaA, alpha, Null.UNIT);
            ijInit(bA, b, Null.UNIT);
            ijInit(C6A, C6, Kelvin.UNIT);
            ijInit(C8A, C8, Kelvin.UNIT);
        }
        else {
            pos = posB;
            for (int i=0; i<q.length; i++) {
                q[i] = sqrtK.toSim(qB[i]);
            }
            ijInit(AB, A, Kelvin.UNIT);
            ijInit(alphaB, alpha, Null.UNIT);
            ijInit(bB, b, Null.UNIT);
            ijInit(C6B, C6, Kelvin.UNIT);
            ijInit(C8B, C8, Kelvin.UNIT);
        }
        ri = space.makeVector();
        rj = space.makeVector();
    }
    
    protected void ijInit(double[] x, double[][] xx, Unit unit) {
        int k = 0;
        for (int i=0; i<4; i++) {
            xx[i][i] = unit.toSim(x[k]);
            k++;
            for (int j=i+1; j<4; j++) {
                xx[i][j] = xx[j][i] = unit.toSim(x[k]);
                k++;
            }
        }
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    public IVector[] gradient(IAtomList atoms) {
        return null;
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return null;
    }
    protected boolean debug = false;

    public double energy(IAtomList atoms) {
        IAtomOriented atom0 = (IAtomOriented)atoms.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)atoms.getAtom(1);
        IVector cm0 = atom0.getPosition();
        IVector cm1 = atom1.getPosition();
        IVector or0 = atom0.getOrientation().getDirection();
        IVector or1 = atom1.getOrientation().getDirection();
        double u = 0;
        boolean checkme = false;
        for (int i=0; i<7; i++) {
            int ii = siteID[i];
            ri.E(cm0);
            ri.PEa1Tv1(pos[i], or0);
            for (int j=0; j<7; j++) {
                int jj = siteID[j];
                rj.E(cm1);
                rj.PEa1Tv1(pos[j], or1);
                double rij2 = rj.Mv1Squared(ri);
                double rij = Math.sqrt(rij2);
                if (rij < 0.9) return Double.POSITIVE_INFINITY;
                if (rij < 1.2) checkme = true;
                double uExp = A[ii][jj]*Math.exp(-alpha[ii][jj]*rij);

                double sum = 1;
                double br = b[ii][jj]*rij;
                double term = 1;
                for (int k=1; k<=6; k++) {
                    term *= br/k;
                    sum += term;
                }
                if (sum==1) return Double.POSITIVE_INFINITY;
                double rij6 = rij2*rij2*rij2;
                double expbr = Math.exp(-br);
                double u6 = -(1-expbr*sum)*C6[ii][jj]/rij6;

                for (int k=7; k<=8; k++) {
                    term *= br/k;
                    sum += term;
                }
                double rij8 = rij6*rij2;
                double u8 = -(1-expbr*sum)*C8[ii][jj]/rij8;

                double uCharge = q[ii]*q[jj]/rij;
                u += uExp + u6 + u8 + uCharge;
                if (debug) {
                    System.out.println(rij+" "+(uExp+u6+u8+uCharge));
                }
                if (Double.isNaN(u)) throw new RuntimeException("oops");
            }
        }
        if (debug) return u;
        if (u<-1000) {
            System.out.println(u);
            debug = true;
            energy(atoms);
            throw new RuntimeException("oops, too much");
        }
        if (checkme && u<10000) {
            System.out.println(u);
        }
        return u;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(IBox box) {
        
    }

    public int nBody() {
        return 2;
    }

    public IVector[][] gradientAndTorque(IAtomList atoms) {
        return null;
    }
    
    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating species = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(species);
        IBox box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(species, 2);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
        ((IAtomOriented)pair.getAtom(0)).getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
        atom1.getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
        P2CO2Hellmann p2 = new P2CO2Hellmann(space, Parameters.B);
        System.out.println("or: "+((IAtomOriented)pair.getAtom(0)).getOrientation().getDirection()+" "+atom1.getOrientation().getDirection());
        for (int i=13; i<=48; i++) {
            atom1.getPosition().setX(0, i*0.25);
            System.out.print(String.format("%5.2f  %20.15e\n", i*0.25, Kelvin.UNIT.fromSim(p2.energy(pair))));
        }
    }
}
