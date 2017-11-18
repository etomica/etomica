
package etomica.models.water;

import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMolecular;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * 4-point potential for water.  Potential parameters are typically defined
 * by a convenience subclass.
 * <p>
 * add a hard core for dielectric virial calculation
 * cannot be used for MC simulation! no nearest boudary condition is used!!!!!
 *
 * @author David Kofke, Andrew Schultz
 */
public class P2WaterTIP4PHardCore extends PotentialMolecular {

    //	private static double A = 600e3; // kcal A^12 / mol
//	private static double C = 610; // kcal A^6 / mol
//	private static double s6 = A/C;
//	private static double s = Math.pow(s6, 1.0/6.0);
//	private static double e = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C/s6*1000))/4.0;
// we can get these parameters from main method
    public double sigma, sigma2;
    protected double epsilon, epsilon4;
    public double hardCore;
    protected Boundary boundary;
    protected final double chargeH;
    protected final double chargeM;
    protected final double chargeMM, chargeMH, chargeHH;
    private Vector P1r, P2r;//dipole position
    private Vector O1P1r, O2P2r;
    private Vector P1P2r;//vector between dipoles
    private Vector mu1Normalized, mu2Normalized;
    private double mu;//magnitude of dipole moment
    private Vector r;
    private Vector dipole1, dipole2;
    private double u_DipoleDipole, uLJDipoleAtCenter;
    private double offset = 0.36794113830914743;//distance between lj site and "P", dipole position
    private double dipoleSeperation = 0.43588227661829493;

    public P2WaterTIP4PHardCore(Space space, double hardCore, double sigma, double epsilon, double chargeM, double chargeH) {
        super(2, space);
        this.hardCore = hardCore;
        this.sigma = sigma;
        sigma2 = sigma * sigma;
        this.epsilon = epsilon;
        epsilon4 = 4 * epsilon;
        chargeMM = chargeM * chargeM;
        chargeMH = chargeM * chargeH;
        chargeHH = chargeH * chargeH;
        this.chargeM = chargeM;
        this.chargeH = chargeH;
        P1r = space.makeVector();
        P2r = space.makeVector();
        O1P1r = space.makeVector();
        O2P2r = space.makeVector();
        P1P2r = space.makeVector();
        mu1Normalized = space.makeVector();
        mu2Normalized = space.makeVector();
        r = space.makeVector();
        dipole1 = space.makeVector();
        dipole2 = space.makeVector();
        mu = -chargeM * dipoleSeperation;//magnitude of dipole moment
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public double energy(IMoleculeList pair) {
        double sum = 0.0;
        IMolecule water1 = pair.getMolecule(0);
        IMolecule water2 = pair.getMolecule(1);

        Vector O1r = water1.getChildList().getAtom(2).getPosition();//H-H-O-M, so O is the third atom
        Vector O2r = water2.getChildList().getAtom(2).getPosition();
        Vector H11r = water1.getChildList().getAtom(0).getPosition();// 1st H on water 1
        Vector H12r = water1.getChildList().getAtom(1).getPosition();// 2nd H on water 1
        Vector H21r = water2.getChildList().getAtom(0).getPosition();// 1st H on water 2
        Vector H22r = water2.getChildList().getAtom(1).getPosition();// 2nd H on water 2
        Vector M1r = water1.getChildList().getAtom(3).getPosition();
        Vector M2r = water2.getChildList().getAtom(3).getPosition();
        r.Ev1Mv2(O1r, O2r);//distance between lj sites

        O1P1r.Ev1Mv2(M1r, O1r);
        O1P1r.normalize();
        mu1Normalized.E(O1P1r);
        O1P1r.TE(offset);
        P1r.Ev1Pv2(O1P1r, O1r);//get position of P1

        O2P2r.Ev1Mv2(M2r, O2r);
        O2P2r.normalize();
        mu2Normalized.E(O2P2r);
        O2P2r.TE(offset);
        P2r.Ev1Pv2(O2P2r, O2r);//get position of P2
        P1P2r.Ev1Mv2(P2r, P1r);//get P1P2r,vector between two dipoles

        double dipoleDistance2 = P1P2r.squared();//correspond to the distance between two dipoles
        P1P2r.normalize();
        //double ulj = sum;
        //boolean printme = r2 > 400;

        // calculate LJ potential at dipole center
        double s2DipoleCenter = sigma2 / dipoleDistance2;
        double s6DipoleCenter = s2DipoleCenter * s2DipoleCenter * s2DipoleCenter;
        uLJDipoleAtCenter = epsilon4 * s6DipoleCenter * (s6DipoleCenter - 1.0);

        double dipoleDistance = Math.sqrt(dipoleDistance2);

        // calculate Dipole-Dipole potential from the other formula
        u_DipoleDipole = mu * mu * (mu1Normalized.dot(mu2Normalized) - 3.0 * mu1Normalized.dot(P1P2r) * P1P2r.dot(mu2Normalized)) / dipoleDistance / dipoleDistance2;

        double r2LJ = r.squared();
        if (r2LJ < hardCore * hardCore) {
            return Double.POSITIVE_INFINITY;
        }

        double s2 = sigma2 / (r2LJ);
        double s6 = s2 * s2 * s2;
        sum += epsilon4 * s6 * (s6 - 1.0);// LJ energy, part of water potential

        double r2 = M1r.Mv1Squared(M2r);
        sum += chargeMM / Math.sqrt(r2);
        r2 = M1r.Mv1Squared(H21r);
        sum += chargeMH / Math.sqrt(r2);
        r2 = M1r.Mv1Squared(H22r);
        sum += chargeMH / Math.sqrt(r2);
        r2 = H11r.Mv1Squared(M2r);
        sum += chargeMH / Math.sqrt(r2);
        r2 = H11r.Mv1Squared(H21r);
        sum += chargeHH / Math.sqrt(r2);
        r2 = H11r.Mv1Squared(H22r);
        sum += chargeHH / Math.sqrt(r2);
        r2 = H12r.Mv1Squared(M2r);
        sum += chargeMH / Math.sqrt(r2);
        r2 = H12r.Mv1Squared(H21r);
        sum += chargeHH / Math.sqrt(r2);
        r2 = H12r.Mv1Squared(H22r);
        sum += chargeHH / Math.sqrt(r2);
//        System.out.println("u_DipoleDipole"+u_DipoleDipole);
//        System.out.println("uLJDipoleAtCenter"+uLJDipoleAtCenter);
//        System.out.println("ulj"+ulj);
//        System.out.println("sum:"+sum);
        return sum;
    }

    public Vector getDipole1() {
        dipole1.E(mu1Normalized);
        dipole1.TE(mu);
        return dipole1;
    }

    public Vector getDipole2() {
        dipole2.E(mu2Normalized);
        dipole2.TE(mu);
        return dipole2;
    }

    public double getULJAtDipoleCenter() {
        return uLJDipoleAtCenter;
    }

    public double getU_DipoleDipole() {
        return u_DipoleDipole;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double getSigma() {
        return sigma;
    }

    public double getEpsilon() {
        return epsilon;
    }
}
