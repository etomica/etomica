package etomica.models.water;

import Jama.Matrix;
import etomica.atom.AtomPair;
import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.math.SpecialFunctions;
import etomica.potential.PotentialN;
import etomica.potential.PotentialPolarizable;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.Electron;
import etomica.units.Kelvin;
import etomica.util.Arrays;

/**
 * GCPM Water potential class.  This class assumes assumes no periodic
 * boundaries exist.  The polarization energy is solved for using matrix
 * inversion rather than iteration, so this class may not be suitable for
 * large systems.
 * 
 * @author Ken
 */
public class PNWaterGCPM extends PotentialN implements PotentialPolarizable {

    public PNWaterGCPM(Space space) {
	    super(space);
	    pair = new AtomPair();
        sigma = 3.69;
        epsilon = Kelvin.UNIT.toSim(110);
        gamma = 12.75;
        chargeH = Electron.UNIT.toSim(0.6113);
        chargeM = Electron.UNIT.toSim(-1.2226);
        core = 4.41; //4.41 = 2.1^2; value according to Cummings
        sigmaM = 0.610;
        sigmaH = 0.455;
        sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        massH = 1.01;
        massO = 16.0;
        totalMass = 18.02;
        sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));
        alphaPol = 1.444;
        
        comWi = space.makeVector();
        comWj = space.makeVector();

        rijVector = space.makeVector();

        work = space.makeVector();
        
        Tunit = space.makeTensor();
        Tij = space.makeTensor();

        Eq = new Matrix[0];
        A = new Matrix[0];
	}   

    public double energy(AtomSet atoms){
        double sum = 0;
        for (int i=0; i<atoms.getAtomCount()-1; i++) {
            pair.atom0 = atoms.getAtom(i);
            for (int j=i+1; j<atoms.getAtomCount(); j++) {
                pair.atom1 = atoms.getAtom(j);
                sum += getNonPolarizationEnergy(pair);
                if (Double.isInfinite(sum)) {
                    return sum;
                }
            }
        }
        
        sum += getPolarizationEnergy(atoms);
        return sum;
    }
    
    /**
     * This returns the pairwise-additive portion of the GCPM potential for a
     * pair of atoms (dispersion + fixed-charge electrostatics)
     */
    public double getNonPolarizationEnergy(AtomSet atoms) {
        AtomWater4P water1 = (AtomWater4P)atoms.getAtom(0);
        AtomWater4P water2 = (AtomWater4P)atoms.getAtom(1);

        IVector O1r = water1.O.getPosition();
        IVector O2r = water2.O.getPosition();
        IVector H11r = water1.H1.getPosition();
        IVector H12r = water1.H2.getPosition();
        IVector H21r = water2.H1.getPosition();
        IVector H22r = water2.H2.getPosition();

        IVector M1r = water1.M.getPosition();
        IVector M2r = water2.M.getPosition();

        double r2 = O1r.Mv1Squared(O2r);
        
        if(r2<=core) {
            return Double.POSITIVE_INFINITY;
        }
        
        double r = Math.sqrt(r2);
        double rOverSigma = r/sigma;
        double sigma2OverR2 = 1/(rOverSigma*rOverSigma);
        double sixOverGamma = 6/gamma;
   
        double sum = epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rOverSigma)) - sigma2OverR2*sigma2OverR2*sigma2OverR2);
        
        r2 = H11r.Mv1Squared(H21r);
        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H22r);
        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H21r);
        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H22r);
        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = M1r.Mv1Squared(H21r);
        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H22r);
        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H11r);
        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H12r);
        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(M2r);
        sum += chargeM*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));

        return sum;
    }

    /**
     * This returns the polarizable portion of the GCPM potential for any
     * number of atoms.
     */
    public double getPolarizationEnergy(AtomSet atoms) {
        
        final int atomCount = atoms.getAtomCount();
        if (Eq.length < atomCount+1) {
            Eq = (Matrix[])Arrays.resizeArray(Eq, atomCount+1);
            A = (Matrix[])Arrays.resizeArray(A, atomCount+1);
        }
        if (Eq[atomCount] == null) {
            Eq[atomCount] = new Matrix(3*atomCount, 1);
            A[atomCount] = new Matrix(3*atomCount, 3*atomCount);
            
            for (int i=0; i<3*atomCount; i++) {
                A[atomCount].set(i, i, 1);
            }
        }
        final Matrix myEq = Eq[atomCount];
        final Matrix myA = A[atomCount];
        for (int i=0; i<3*atomCount; i++) {
            myEq.set(i, 0, 0);
        }
        
        /*
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         */

        for (int i=0; i<atoms.getAtomCount(); i++) {
            AtomWater4P atomi = (AtomWater4P)atoms.getAtom(i);
            IVector O1r = atomi.O.getPosition();
            IVector H11r = atomi.H1.getPosition();
            IVector H12r = atomi.H2.getPosition();

            comWi.Ea1Tv1(massH, H11r);
            comWi.PEa1Tv1(massO, O1r);
            comWi.PEa1Tv1(massH, H12r);
            comWi.TE(1.0/totalMass);
            
            for (int j=0; j<atoms.getAtomCount(); j++) {
                if  (i == j) continue;
                AtomWater4P atomj = (AtomWater4P)atoms.getAtom(j);
            
                double comWtoH1 = Math.sqrt(comWi.Mv1Squared(atomj.H1.getPosition()));
                double comWtoH2 = Math.sqrt(comWi.Mv1Squared(atomj.H2.getPosition()));
                double comWtoM = Math.sqrt(comWi.Mv1Squared(atomj.M.getPosition()));
    
                double fac = chargeH/(comWtoH1*comWtoH1*comWtoH1)*((1-SpecialFunctions.erfc(comWtoH1/sqrtHMsigmas))
                        -Math.sqrt(2)*comWtoH1/sqrtPiHMsigmas*Math.exp(-comWtoH1*comWtoH1/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
                work.Ev1Mv2(comWi, atomj.H1.getPosition());
                work.TE(fac);
                myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+work.x(0));
                myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+work.x(1));
                myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+work.x(2));
    
                fac = chargeH/(comWtoH2*comWtoH2*comWtoH2)*((1-SpecialFunctions.erfc(comWtoH2/sqrtHMsigmas))
                        -Math.sqrt(2)*comWtoH2/sqrtPiHMsigmas*Math.exp(-comWtoH2*comWtoH2/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
                work.Ev1Mv2(comWi, atomj.H2.getPosition());
                work.TE(fac);
                myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+work.x(0));
                myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+work.x(1));
                myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+work.x(2));
    
                fac = chargeM/(comWtoM*comWtoM*comWtoM)*((1-SpecialFunctions.erfc(comWtoM/(2*sigmaM)))
                        -Math.sqrt(2)*comWtoM/sqrtPiMMsigmas*Math.exp(-comWtoM*comWtoM/(4*sigmaM*sigmaM)));
                work.Ev1Mv2(comWi, atomj.M.getPosition());
                work.TE(fac);
                myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+work.x(0));
                myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+work.x(1));
                myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+work.x(2));
                
                if (i<j) {
                    comWj.Ea1Tv1(massH, atomj.H1.getPosition());
                    comWj.PEa1Tv1(massO, atomj.O.getPosition());
                    comWj.PEa1Tv1(massH, atomj.H2.getPosition());
                    comWj.TE(1.0/totalMass);
                    
                    rijVector.Ev1Mv2(comWj,comWi);  // is this the correct direction? kmb, 8/7/06
                    
                    double r12 = Math.sqrt(rijVector.squared());


                    double f = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)) + (r12*r12*r12)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
                    
                    double g = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
                    
                    // Filling the unit matrix I
                    Tij.Ev1v2(rijVector,rijVector);
                    
                    Tij.TE(3*f/(r12*r12));
                    
                    Tunit.E(g);
                    
                    Tij.ME(Tunit);
                    Tij.TE(1/(r12*r12*r12));
                    
                    //Try matrix inversion solution with Jama library
                            
                    Tij.TE(alphaPol);
                    
                    int mOffset = i*3;
                    int nOffset = j*3;
                    for (int m=0; m<3; m++) {
                        for (int n=0; n<3; n++) {
                            myA.set(mOffset+m, nOffset+n, -Tij.component(m, n));
                            myA.set(nOffset+n, mOffset+m, -Tij.component(n, m));
                        }
                    }
                }
            }
            
        }
        
        //x here represents P (almost).
        //For x to be P, the A of the Ax=b actually needs an extra factor of
        //alphaPol.  We'll add that bit in when we calculate UpolAtkins.  

        Matrix x = myA.solve(myEq);

        if (false) {
            // this is (mathematically) what we want.  But Jama is slow.
            UpolAtkins = -0.5*(x.transpose().times(myEq)).get(0,0)*alphaPol;
        }
        else {
            UpolAtkins = 0;
            for (int i=0; i<3*atomCount; i++) {
                UpolAtkins += x.get(i,0)*myEq.get(i,0);
            }
            UpolAtkins *= -0.5*alphaPol;
        }

        // only needed for more complicated Eq8 from Cummings paper 
        if (false) {
            
            // for the sake of clarity (over perf), just multiply x by alphaPol
            // (see comment above about A lacking alphaPol)
            x.timesEquals(alphaPol);
            Matrix Ep = myA.times(x).minus(x);
            Ep.timesEquals(-1/alphaPol);

            double x2NormF = x.normF();
            double UpolEquation8 = 2*UpolAtkins -0.5*(x.transpose().times(Ep).get(0,0))+(0.5/alphaPol)*(x2NormF*x2NormF);

            if (Math.abs(UpolAtkins-UpolEquation8) > 1.e-6) {
                throw new RuntimeException("oops "+UpolAtkins+" "+UpolEquation8);
            }
        }
        
        return UpolAtkins;
    }
    
    public double getLastPolarizationEnergy() {
        return UpolAtkins;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public void setBox(Box box) {
    }

    private static final long serialVersionUID = 1L;
    protected final AtomPair pair;
    private final double sigma;
    private final double epsilon, gamma;
    private final double chargeH, chargeM;
    private final double core; // = 4.41; //4.41 = 2.1^2; value according to Cummings
    private Matrix[] Eq, A;
    private IVector comWi, comWj;
    private final IVector rijVector;
    private final IVector work;
    private final Tensor Tunit, Tij;
    private final double sigmaM;
    private final double sigmaH;
    private final double sqrtHMsigmas;
    private final double massH;
    private final double massO;
    private final double totalMass;
    private final double sqrtPiHMsigmas;
    private final double sqrtPiMMsigmas;
    private final double alphaPol;
    private double UpolAtkins;
}