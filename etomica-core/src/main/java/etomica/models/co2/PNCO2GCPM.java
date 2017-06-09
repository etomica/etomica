/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.co2;

import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.IElement;
import etomica.chem.elements.Oxygen;
import etomica.config.IConformation;
import etomica.models.co2.P2CO2Hellmann.Parameters;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.potential.IPotentialMolecular;
import etomica.potential.P3AxilrodTeller;
import etomica.potential.PotentialMolecular;
import etomica.potential.PotentialPolarizable;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresHetero;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Electron;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;

/**
 * GCPM CO2 potential class (GCPCDO).  This potential was described in
 * 
 * http://dx.doi.org/10.1063/1.3519022
 * 
 * and checked against Persson's fortran implementation
 * 
 * @author Ken, Andrew and Dave
 */
public class PNCO2GCPM extends PotentialMolecular implements PotentialPolarizable {

    public PNCO2GCPM(Space space) {
        this(space, Integer.MAX_VALUE);
    }
    
    public PNCO2GCPM(Space space, int nBody) {
	    super(nBody, space);
	    pair = new MoleculePair();
        sigmaC = 3.193;
        sigmaO = sigmaC*1.0483; //paper says 3.347;
        coreFac = 0.57*0.57;
        epsilonC = Kelvin.UNIT.toSim(71.34);
        epsilonO = Kelvin.UNIT.toSim(67.72);
        tauO = 0.6100;
        tauC = tauO/1.0483;  // paper reports 0.5819, fortran uses this instead
        sqrtCOtau = Math.sqrt(2*(tauC*tauC+tauO*tauO));
        double[] t = new double[]{tauC,tauO,tauO};
        double[] s = new double[]{sigmaC,sigmaO,sigmaO};
        double[] e = new double[]{epsilonC,epsilonO,epsilonO};
        tauAll = new double[3][3];
        sigmaAll = new double[3][3];
        epsilonAll = new double[3][3];
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                sigmaAll[i][j] = 0.5*(s[i]+s[j]);
                epsilonAll[i][j] = 2*e[i]*e[j]/(e[i]+e[j]);
                tauAll[i][j] = Math.sqrt(2*(t[i]*t[i]+t[j]*t[j]));
            }
        }
        gamma = 15.5;
        chargeO = Electron.UNIT.toSim(-0.3321);
        chargeC = Electron.UNIT.toSim(+0.6642);
        chargeAll = new double[]{chargeC,chargeO,chargeO};
        sqrtPiCOtau = Math.sqrt(Math.PI*(tauC*tauC+tauO*tauO));
        sqrtPiCCtau = Math.sqrt(Math.PI*(2*tauC*tauC));
        alphaPar = 4.05;
        alphaPerp = 1.95;
        
        oldMu = space.makeVector();
        shift = space.makeVector();

        rijVector = space.makeVector();

        work = space.makeVector();
        
        Tunit = space.makeTensor();
        Tunit.E(new double[][]{{1,0,0},{0,1,0},{0,0,1}});
        Tij = space.makeTensor();

        Eq = new Vector[0];
        Ep = new Vector[0];
        mu = new Vector[0];
        component = Component.FULL;
	}
    
    public void setComponent(Component comp) {
        component = comp;
    }

    protected boolean oops = false;
    public double energy(IMoleculeList atoms){
        double sum = 0;
        if (component != Component.INDUCTION) {
            for (int i=0; i<atoms.getMoleculeCount()-1; i++) {
                pair.atom0 = atoms.getMolecule(i);
                for (int j=i+1; j<atoms.getMoleculeCount(); j++) {
                    pair.atom1 = atoms.getMolecule(j);
                    sum += getNonPolarizationEnergy(pair);
                    if (Double.isInfinite(sum)) {
                        return sum;
                    }
                }
            }
        }
        if (component != Component.TWO_BODY) {
            sum += getPolarizationEnergy(atoms);
        }
        if (!oops && Double.isNaN(sum)) {
            oops = true;
            energy(atoms);
            throw new RuntimeException("oops NaN");
        }
        return sum;
    }
    
    public static boolean debugme = false;
    
    /**
     * This returns the pairwise-additive portion of the GCPM potential for a
     * pair of atoms (dispersion + fixed-charge electrostatics)
     */
    public double getNonPolarizationEnergy(IMoleculeList molecules) {
        IAtomList water1Atoms = molecules.getMolecule(0).getChildList();
        IAtomList water2Atoms = molecules.getMolecule(1).getChildList();

        Vector C1r = water1Atoms.getAtom(0).getPosition();
        Vector C2r = water2Atoms.getAtom(0).getPosition();
        
        work.Ev1Mv2(C1r, C2r);
        shift.Ea1Tv1(-1,work);
		boundary.nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;
        
        double r2 = work.squared();
        
        if(r2<=sigmaC*coreFac) {
            return Double.POSITIVE_INFINITY;
        }

        Vector O11r = water1Atoms.getAtom(1).getPosition();
        Vector O12r = water1Atoms.getAtom(2).getPosition();
        Vector O21r = water2Atoms.getAtom(1).getPosition();
        Vector O22r = water2Atoms.getAtom(2).getPosition();

        double sum =0;
        if (zeroShift) {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    r2 = water1Atoms.getAtom(i).getPosition().Mv1Squared(water2Atoms.getAtom(j).getPosition());
                    double r = Math.sqrt(r2);
                    double rOverSigma = r/sigmaAll[i][j];
                    double sigma2OverR2 = 1/(rOverSigma*rOverSigma);
                    if (1/sigma2OverR2 < coreFac) return Double.POSITIVE_INFINITY;
                    double sixOverGamma = 6/gamma;
                    sum += epsilonAll[i][j]/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rOverSigma)) - sigma2OverR2*sigma2OverR2*sigma2OverR2);//exp-6 potential(Udisp)
                }
            }
        }
        else {
            for (int i=0; i<3; i++) {
                for (int j=0; j<3; j++) {
                    Vector r1 = water1Atoms.getAtom(i).getPosition();
                    shift.PE(r1);
                    r2 = water2Atoms.getAtom(j).getPosition().Mv1Squared(shift);
                    shift.ME(r1);
                    
                    double r = Math.sqrt(r2);
                    double rOverSigma = r/sigmaAll[i][j];
                    double sigma2OverR2 = 1/(rOverSigma*rOverSigma);
                    if (1/sigma2OverR2 < coreFac) return Double.POSITIVE_INFINITY;
                    double sixOverGamma = 6/gamma;
               
                    sum += epsilonAll[i][j]/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rOverSigma)) - sigma2OverR2*sigma2OverR2*sigma2OverR2);//exp-6 potential(Udisp)
                }
            }
        }

        if (zeroShift){
	        r2 = O11r.Mv1Squared(O21r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));

	        r2 = O11r.Mv1Squared(O22r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	
	        r2 = O12r.Mv1Squared(O21r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	
	        r2 = O12r.Mv1Squared(O22r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	        
	        r2 = C1r.Mv1Squared(O21r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        r2 = C1r.Mv1Squared(O22r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        r2 = C2r.Mv1Squared(O11r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        r2 = C2r.Mv1Squared(O12r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));

	        r2 = C1r.Mv1Squared(C2r);
	        sum += chargeC*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauC)));
        }
        else {
        	shift.PE(O11r);
        	r2 = O21r.Mv1Squared(shift);
        	shift.ME(O11r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	        
	        shift.PE(O11r);
        	r2 = O22r.Mv1Squared(shift);
        	shift.ME(O11r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	
	        shift.PE(O12r);
        	r2 = O21r.Mv1Squared(shift);
        	shift.ME(O12r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	
	        shift.PE(O12r);
        	r2 = O22r.Mv1Squared(shift);
        	shift.ME(O12r);
	        sum += chargeO*chargeO/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauO)));
	        
	        shift.PE(C1r);
        	r2 = O21r.Mv1Squared(shift);
        	shift.ME(C1r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        shift.PE(C1r);
        	r2 = O22r.Mv1Squared(shift);
        	shift.ME(C1r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        shift.PE(O11r);
        	r2 = C2r.Mv1Squared(shift);
        	shift.ME(O11r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        shift.PE(O12r);
        	r2 = C2r.Mv1Squared(shift);
        	shift.ME(O12r);
	        sum += chargeO*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/sqrtCOtau));
	
	        shift.PE(C1r);
        	r2 = C2r.Mv1Squared(shift);
        	shift.ME(C1r);
	        sum += chargeC*chargeC/Math.sqrt(r2)*(1-org.apache.commons.math3.special.Erf.erfc(Math.sqrt(r2)/(2*tauC)));
        	
        }
        return sum;
    }

    /**
     * This returns the polarizable portion of the GCPM potential for any
     * number of atoms.
     */
    public double getPolarizationEnergy(IMoleculeList molecules) {
        
        final int atomCount = molecules.getMoleculeCount();
        if (Eq.length < atomCount+1) {
            int oldSize = Eq.length;
            Eq = (Vector[])etomica.util.Arrays.resizeArray(Eq, atomCount);
            Ep = (Vector[])etomica.util.Arrays.resizeArray(Ep, atomCount);
            mu= (Vector[])etomica.util.Arrays.resizeArray(mu, atomCount);
            for (int i=oldSize; i<atomCount; i++) {
                Eq[i] = space.makeVector();
                Ep[i] = space.makeVector();
                mu[i] = space.makeVector();
            }
        }
        for (int i=0; i<atomCount; i++) {
            Eq[i].E(0);
            mu[i].E(0);
            Ep[i].E(0);
        }
        
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList iLeafAtoms = molecules.getMolecule(i).getChildList();
            Vector C1r = iLeafAtoms.getAtom(0).getPosition();

            for (int j=0; j<molecules.getMoleculeCount(); j++) {
                if  (i == j) continue;
                IAtomList jLeafAtoms = molecules.getMolecule(j).getChildList();
                Vector Cjr = jLeafAtoms.getAtom(0).getPosition();
                Vector Oj1r = jLeafAtoms.getAtom(1).getPosition();
                Vector Oj2r = jLeafAtoms.getAtom(2).getPosition();
                
                work.Ev1Mv2(C1r, Cjr);
                shift.Ea1Tv1(-1,work);
                boundary.nearestImage(work);
                shift.PE(work);
                final boolean zeroShift = shift.squared() < 0.1;
                double rCtoO1,rCtoO2,rCtoC;
                
                if (zeroShift){
                
                    rCtoO1 = Math.sqrt(C1r.Mv1Squared(Oj1r));
                    rCtoO2 = Math.sqrt(C1r.Mv1Squared(Oj2r));
                    rCtoC = Math.sqrt(C1r.Mv1Squared(Cjr));
                }
                else {
                    shift.PE(C1r);
                    rCtoO1 = Math.sqrt(Oj1r.Mv1Squared(shift));
                    shift.ME(C1r);
                    
                    shift.PE(C1r);
                    rCtoO2 = Math.sqrt(Oj2r.Mv1Squared(shift));
                    shift.ME(C1r);
                    
                    shift.PE(C1r);
                    rCtoC = Math.sqrt(Cjr.Mv1Squared(shift));
                    shift.ME(C1r);
                }
                
                
                
                // For molecules that are far apart, fac=chargeX/comWtoX^3, but we add up
                // facs for H and M, which mostly cancel each other out, so we lose quite 
                // a bit of precision (~2-3 digits).
                double fac = chargeO/(rCtoO1*rCtoO1*rCtoO1)*((1-org.apache.commons.math3.special.Erf.erfc(rCtoO1/sqrtCOtau))
                        -Math.sqrt(2)*rCtoO1/sqrtPiCOtau*Math.exp(-rCtoO1*rCtoO1/(2*(tauC*tauC+tauO*tauO))));
                work.Ev1Mv2(C1r, Oj1r);
                work.PE(shift);
                work.TE(fac);
                Eq[i].PE(work);
    
                fac = chargeO/(rCtoO2*rCtoO2*rCtoO2)*((1-org.apache.commons.math3.special.Erf.erfc(rCtoO2/sqrtCOtau))
                        -Math.sqrt(2)*rCtoO2/sqrtPiCOtau*Math.exp(-rCtoO2*rCtoO2/(2*(tauC*tauC+tauO*tauO))));
                work.Ev1Mv2(C1r, Oj2r);
                work.PE(shift);
                work.TE(fac);
                Eq[i].PE(work);
    
                fac = chargeC/(rCtoC*rCtoC*rCtoC)*((1-org.apache.commons.math3.special.Erf.erfc(rCtoC/(2*tauC)))
                        -Math.sqrt(2)*rCtoC/sqrtPiCCtau*Math.exp(-rCtoC*rCtoC/(4*tauC*tauC)));
                work.Ev1Mv2(C1r, Cjr);
                work.PE(shift);
                work.TE(fac);
                Eq[i].PE(work);
            }
        }

        int maxIter = 550;
        double mixIter = 0.9;
        double sqrtpi = Math.sqrt(Math.PI);
for (int iter=0; iter<maxIter; iter++) {
        double sumDeltaMu = 0;
        double sumMu = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList iLeafAtoms = molecules.getMolecule(i).getChildList();
            Vector C1r = iLeafAtoms.getAtom(0).getPosition();

            work.Ev1Mv2(C1r,iLeafAtoms.getAtom(1).getPosition());
            work.normalize();
            Ep[i].PE(Eq[i]);
            double cosTheta = Math.abs(work.dot(Ep[i])/Math.sqrt(Ep[i].squared()));
            oldMu.E(mu[i]);
            mu[i].Ea1Tv1(alphaPerp + cosTheta*(alphaPar-alphaPerp), Ep[i]);
            mu[i].TE(mixIter);
            mu[i].PEa1Tv1(1-mixIter, oldMu);
            sumDeltaMu += mu[i].Mv1Squared(oldMu);
            sumMu += mu[i].squared();
        }
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            Ep[i].E(0);
        }

        double tauK = tauAll[0][0];
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            IAtomList iLeafAtoms = molecules.getMolecule(i).getChildList();
            Vector C1r = iLeafAtoms.getAtom(0).getPosition();

            for (int j=i+1; j<molecules.getMoleculeCount(); j++) {
                IAtomList jLeafAtoms = molecules.getMolecule(j).getChildList();
                Vector Cjr = jLeafAtoms.getAtom(0).getPosition();
                
                rijVector.Ev1Mv2(C1r, Cjr);
        		boundary.nearestImage(rijVector);
        		
        		double r2 = rijVector.squared();
                     
                if (r2 < coreFac*sigmaC) {
                    UpolAtkins = Double.NaN;
                    return UpolAtkins;
                }
                
                double r1 = Math.sqrt(r2);
                
                double erf = (1-org.apache.commons.math3.special.Erf.erfc(r1/(tauK)));
                double exp = Math.exp(-r2/(tauK*tauK));
                
                double prefac = (r1/(0.5*tauK*sqrtpi))*exp;
                
                double postfac = prefac * 0.666666666666666666666 * r2 / (tauK*tauK);
                
                double fr = erf - prefac;

                double fpr = fr - postfac;
                
                Ep[i].PEa1Tv1(-fr/(r1*r2), mu[j]);
                
                Ep[i].PEa1Tv1(3*rijVector.dot(mu[j])*fpr/(r2*r2*r1), rijVector);
                
                Ep[j].PEa1Tv1(-fr/(r1*r2), mu[i]);
                Ep[j].PEa1Tv1(3*rijVector.dot(mu[i])*fpr/(r2*r2*r1), rijVector);
            }
        }
        
        if (debugme) {
            for (int i=0; i<molecules.getMoleculeCount(); i++) {
                System.out.println(iter+" "+i+" "+Ep[i]+" "+mu[i]);
            }
        }

        if (sumDeltaMu < 1e-20) break;
        if (iter==maxIter-1) {
            System.err.println("we were unable to converge");
            System.err.println("sumDeltaMu "+sumDeltaMu);
            System.err.println("sumMu "+sumMu);
        }
}
        UpolAtkins = 0;
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            UpolAtkins += Eq[i].dot(mu[i]);
        }
        UpolAtkins *= -0.5;
        if (!debugme && Double.isNaN(UpolAtkins)) {
            debugme = true;
            getPolarizationEnergy(molecules);
            throw new RuntimeException("oops");
        }
        //x here represents P (almost).
        //For x to be P, the A of the Ax=b actually needs an extra factor of
        //alphaPol.  We'll add that bit in when we calculate UpolAtkins.  
        return UpolAtkins;
    }
    
    public double getLastPolarizationEnergy() {
        return UpolAtkins;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public void setBox(Box box) {
    	boundary = box.getBoundary();
    }
    
    public P3GCPMAxilrodTeller makeAxilrodTeller() {
        return new P3GCPMAxilrodTeller(space);
    }

    protected final MoleculePair pair;
    protected Boundary boundary;
    protected final double sigmaC, sigmaO;
    protected final double[][] sigmaAll, epsilonAll;
    protected final double epsilonC, epsilonO, gamma;
    protected final double chargeO, chargeC;
    protected final double[] chargeAll;
    protected final double tauC, tauO; // sigma for water
    protected final double[][] tauAll;
    protected final double coreFac;
    protected Vector[] Eq, Ep, mu;
    protected Vector oldMu;
    protected final Vector rijVector;
    protected final Vector work, shift;
    protected final Tensor Tunit, Tij;
    protected final double sqrtCOtau;
    protected final double sqrtPiCOtau;
    protected final double sqrtPiCCtau;
    protected final double alphaPerp, alphaPar;
    protected double UpolAtkins;
    protected Component component;
    
    public enum Component { TWO_BODY, INDUCTION, FULL }
    
    public class P3GCPMAxilrodTeller implements IPotentialMolecular {

        protected final Vector rij, rik, rjk;
        protected final Vector bveci, bvecj, bveck;
        protected final double[] cosg;
        protected final Vector norm;
        protected final Vector xveci, xvecj, xveck;
        protected final Vector yveci, yvecj, yveck;
        protected final double[] xx, yy, zz;
        public static final double dpolx = 1.95, anx = 2.1;
        public final double nufac0;
        
        public P3GCPMAxilrodTeller(Space space) {
            rij = space.makeVector();
            rik = space.makeVector();
            rjk = space.makeVector();
            cosg = new double[4];
            norm = space.makeVector();
            bveci = space.makeVector();
            bvecj = space.makeVector();
            bveck = space.makeVector();
            xveci = space.makeVector();
            xvecj = space.makeVector();
            xveck = space.makeVector();
            yveci = space.makeVector();
            yvecj = space.makeVector();
            yveck = space.makeVector();
            xx = new double[3];
            yy = new double[3];
            zz = new double[3];
            nufac0 = Kelvin.UNIT.toSim(2.52e4);
        }
        
        public double getRange() {
            return Double.POSITIVE_INFINITY;
        }

        public void setBox(Box box) {
            
        }

        public int nBody() {
            return 3;
        }

        public double energy(IMoleculeList molecules) {
            double bfac = 1.0/(2*1.161);
            IAtomList atomsi = molecules.getMolecule(0).getChildList();
            IAtomList atomsj = molecules.getMolecule(1).getChildList();
            IAtomList atomsk = molecules.getMolecule(2).getChildList();
            Vector ri = atomsi.getAtom(0).getPosition();
            Vector rj = atomsj.getAtom(0).getPosition();
            Vector rk = atomsk.getAtom(0).getPosition();
            rij.Ev1Mv2(rj, ri);
            rik.Ev1Mv2(rk, ri);
            rjk.Ev1Mv2(rk, rj);
            double drij2 = rij.squared();
            double drik2 = rik.squared();
            double drjk2 = rjk.squared();
            double drij = Math.sqrt(drij2);
            double drik = Math.sqrt(drik2);
            double drjk = Math.sqrt(drjk2);
            double drij3 = drij2*drij;
            double drik3 = drik2*drik;
            double drjk3 = drjk2*drjk;
            rij.TE(1/drij);
            rik.TE(1/drik);
            rjk.TE(1/drjk);
            cosg[0] = -rij.dot(rik);
            cosg[1] = -rij.dot(rjk);
            cosg[2] = rjk.dot(rik);
            cosg[3] = cosg[0];
            bveci.Ev1Mv2(atomsi.getAtom(1).getPosition(), atomsi.getAtom(2).getPosition());
            bveci.TE(bfac);
            bvecj.Ev1Mv2(atomsj.getAtom(1).getPosition(), atomsj.getAtom(2).getPosition());
            bvecj.TE(bfac);
            bveck.Ev1Mv2(atomsk.getAtom(1).getPosition(), atomsk.getAtom(2).getPosition());
            bveck.TE(bfac);
            norm.E(rij);
            norm.XE(rik);
            norm.normalize();

            xveci.Ev1Pv2(rij, rik);
            xveci.normalize();
            yveci.E(norm);
            yveci.XE(xveci);

            xvecj.Ev1Mv2(rjk,rij);
            xvecj.normalize();
            yvecj.E(norm);
            yvecj.XE(xvecj);
            
            xveck.Ev1Pv2(rjk, rik);
            xveck.TE(-1);
            yveck.E(norm);
            yveck.XE(xveck);

            double eadd = 1;
            for (int a=0; a<3; a++) {
                xx[a] = Math.sqrt((1 + cosg[a])*(1 + cosg[a+1])) +
                        Math.sqrt((1 - cosg[a])*(1 - cosg[a+1])) * 0.5;
                eadd *= xx[a];
            }
            double apol = anx * Math.abs(bveci.dot(xveci));
            double polix = dpolx + apol;
            
            apol = anx * Math.abs(bvecj.dot(xvecj));
            double poljx = dpolx + apol;
            
            apol = anx * Math.abs(bveck.dot(xveck));
            double polkx = dpolx + apol;

            double nufac   = polkx*poljx*polix;
            eadd *= nufac;
            double u = eadd / (drij3*drik3*drjk3);

            // (z, z, z) matrix element
            // the polarizability is here the projection on the normal to the 
            // intermolecular plane. 
            apol = anx * Math.abs(bveck.dot(norm));
            double poliz = dpolx + apol;
            
            apol = anx * Math.abs(bvecj.dot(norm));
            double poljz = dpolx + apol;
            
            apol = anx * Math.abs(bveci.dot(norm));
            double polkz = dpolx + apol;

            nufac   = polkz*poljz*poliz;

            u += nufac / (drij3*drik3*drjk3);

            // (y, y, y) matrix element
            // the polarizability is here the projection on the axis orthogonal both
            // to the bisector axis, and to the normal of the intermolecular plane.
            eadd = 1;
            for (int a=0; a<3; a++) {
                yy[a] = -Math.sqrt((1 - cosg[a])*(1 - cosg[a+1])) -
                        Math.sqrt((1 + cosg[a])*(1 + cosg[a+1])) * 0.5;
                eadd *= yy[a];
            }

            apol = anx * Math.abs(bveci.dot(yveci));
            double poliy = dpolx + apol;
            
            apol = anx * Math.abs(bvecj.dot(yvecj));
            double poljy = dpolx + apol;
            
            apol = anx * Math.abs(bveck.dot(yveck));
            double polky = dpolx + apol;

            nufac   = polkx*poljx*polix;
            u += eadd*nufac / (drij3*drik3*drjk3);

            // here come the mixed matrix elements. six in total. three double x's and
            // three double y's.
            // (x, x, y)
            double eprd = (Math.sqrt((1 + cosg[1])*(1 - cosg[2])) - 
                            Math.sqrt((1 - cosg[1])*(1 + cosg[2])) * 0.5) * 
                           (-Math.sqrt((1 + cosg[0])*(1 - cosg[2])) + 
                            Math.sqrt((1 - cosg[0])*(1 + cosg[2])) * 0.5);
            nufac= polix*poljx*polky;
            eadd = xx[0]*eprd*nufac;

            // (x, y, x)
            eprd = (Math.sqrt((1 + cosg[0])*(1 - cosg[1])) - 
                    Math.sqrt((1 - cosg[0])*(1 + cosg[1])) * 0.5) * 
                  (-Math.sqrt((1 + cosg[2])*(1 - cosg[1])) + 
                    Math.sqrt((1 - cosg[2])*(1 + cosg[1])) * 0.5);
            nufac= polix*poljy*polkx;
            eadd += xx[2]*eprd*nufac;

            // (y, x, x)
            eprd = (Math.sqrt((1 + cosg[2])*(1 - cosg[0])) -
                    Math.sqrt((1 - cosg[2])*(1 + cosg[0])) * 0.5) *
                  (-Math.sqrt((1 + cosg[1])*(1 - cosg[0])) +
                    Math.sqrt((1 - cosg[1])*(1 + cosg[0])) * 0.5);
            nufac= poliy*poljx*polkx;
            eadd += xx[1]*eprd*nufac;
            u  += eadd / (drij3*drik3*drjk3);

            // the double y's.
            // (y, y, x)
            eprd = (Math.sqrt((1.0 + cosg[2])*(1 - cosg[0])) -
                    Math.sqrt((1.0 - cosg[2])*(1 + cosg[0])) * 0.5) *
                  (-Math.sqrt((1.0 + cosg[2])*(1 - cosg[1])) +
                    Math.sqrt((1.0 - cosg[2])*(1 + cosg[1])) * 0.5);
            nufac= poliy*poljy*polkx;
            eadd = yy[0]*eprd*nufac;

            // (y, x, y)
            eprd = (Math.sqrt((1.0 + cosg[1])*(1 - cosg[2])) -
                    Math.sqrt((1.0 - cosg[1])*(1 + cosg[2])) * 0.5) *
                  (-Math.sqrt((1.0 + cosg[1])*(1 - cosg[0])) +
                    Math.sqrt((1.0 - cosg[1])*(1 + cosg[0])) * 0.5);
            nufac= poliy*poljx*polky;
            eadd += yy[2]*eprd*nufac;

            // (x, y, y)
            eprd = (Math.sqrt((1.0 + cosg[0])*(1 - cosg[1])) -
                    Math.sqrt((1.0 - cosg[0])*(1 + cosg[1])) * 0.5) *
                  (-Math.sqrt((1.0 + cosg[0])*(1 - cosg[2])) +
                    Math.sqrt((1.0 - cosg[0])*(1 + cosg[2])) * 0.5);
            nufac= polix*poljy*polky;
            eadd = eadd + yy[1]*eprd*nufac;

            u += eadd / (drij3*drik3*drjk3);
            
            if (Double.isNaN(u)) {
                energy(molecules);
                throw new RuntimeException("oops "+u);
            }
            return u*nufac0;
        }
    }
    
    public static void main2(String[] args) {
        double x = 0;
        double z = 4.;
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesSpheresHetero speciesCO2 = new SpeciesSpheresHetero(space, new IElement[]{Carbon.INSTANCE, Oxygen.INSTANCE});
        speciesCO2.setChildCount(new int[]{1,2});
        speciesCO2.setConformation(new IConformation() {
            
            public void initializePositions(IAtomList atomList) {
                atomList.getAtom(0).getPosition().E(0);
                atomList.getAtom(1).getPosition().setX(0,1.161);
                atomList.getAtom(2).getPosition().setX(0,-1.161);
            }
        });
        sim.addSpecies(speciesCO2);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 2);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IMolecule mol0 = box.getMoleculeList().getMolecule(0);
        IMolecule mol1 = box.getMoleculeList().getMolecule(1);
        
        mol0.getChildList().getAtom(0).getPosition().E(space.makeVector(new double[]{0.000000,0,0.000000 }));
        mol0.getChildList().getAtom(1).getPosition().E(space.makeVector(new double[]{-1.161,0,0 }));
        mol0.getChildList().getAtom(2).getPosition().E(space.makeVector(new double[]{1.161,0,0 }));
        mol1.getChildList().getAtom(0).getPosition().E(space.makeVector(new double[]{ 0,0,z }));
        mol1.getChildList().getAtom(1).getPosition().E(space.makeVector(new double[]{ -1.161,0,z}));
        mol1.getChildList().getAtom(2).getPosition().E(space.makeVector(new double[]{ 1.161,0,z }));
        
//        space.makeVector(new double[]{ 1.000000,-11.000000,-5.000000 }) 
//        space.makeVector(new double[]{ 0.732908,-10.699688,-3.910782 }) 
//        space.makeVector(new double[]{ 1.267092,-11.300312,-6.089218 }) 
        
//        MoleculeActionTranslateTo translator = new MoleculeActionTranslateTo(space);
//        translator.setDestination(space.makeVector(new double[]{x,0,z}));
//        translator.actionPerformed(mol1);
        PNCO2GCPM p2 = new PNCO2GCPM(space);
        p2.setBox(box);
        IMoleculeList molecules = box.getMoleculeList();
        double u = p2.energy(molecules);
        System.out.println(u);
        
        sim = new Simulation(space);
        SpeciesSpheresRotating species2CO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(species2CO2);
        box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(species2CO2, 2);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
        atom1.getPosition().E(space.makeVector(new double[]{x,0,z}));
//        ((IAtomOriented)pair.getAtom(0)).getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
//        IVectorMutable o1 = space.makeVector(new double[]{-1,0,0});
//        atom1.getOrientation().setDirection(o1);
        P2CO2Hellmann p2H = new P2CO2Hellmann(space, Parameters.B);
        double uH = p2H.energy(pair);
        System.out.println("Hellmann: "+uH);
        
    }
    
    public static void main(String[] args) {
        double nufac0 = Kelvin.UNIT.toSim(2.52e4);
        double nufac = 9.0/16.0*ElectronVolt.UNIT.toSim(13.7);
        System.out.println(nufac0+" "+nufac);
        double x1 = 6;
        double z1 = 5.;
        double y2 = 7.;
        double z2 = -2;
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesSpheresHetero speciesCO2 = new SpeciesSpheresHetero(space, new IElement[]{Carbon.INSTANCE, Oxygen.INSTANCE});
        speciesCO2.setChildCount(new int[]{1,2});
        speciesCO2.setConformation(new IConformation() {
            
            public void initializePositions(IAtomList atomList) {
                atomList.getAtom(0).getPosition().E(0);
                atomList.getAtom(1).getPosition().setX(0,1.161);
                atomList.getAtom(2).getPosition().setX(0,-1.161);
            }
        });
        sim.addSpecies(speciesCO2);
        Box box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(speciesCO2, 3);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IMolecule mol0 = box.getMoleculeList().getMolecule(0);
        IMolecule mol1 = box.getMoleculeList().getMolecule(1);
        IMolecule mol2 = box.getMoleculeList().getMolecule(2);
        
        mol0.getChildList().getAtom(0).getPosition().E(space.makeVector(new double[]{0.000000,0,0.000000 }));
        mol0.getChildList().getAtom(1).getPosition().E(space.makeVector(new double[]{-1.161,0,0 }));
        mol0.getChildList().getAtom(2).getPosition().E(space.makeVector(new double[]{1.161,0,0 }));
        mol1.getChildList().getAtom(0).getPosition().E(space.makeVector(new double[]{ x1,0,z1 }));
        mol1.getChildList().getAtom(1).getPosition().E(space.makeVector(new double[]{ x1,0,z1-1.161}));
        mol1.getChildList().getAtom(2).getPosition().E(space.makeVector(new double[]{ x1,0,z1+1.161 }));
        mol2.getChildList().getAtom(0).getPosition().E(space.makeVector(new double[]{ 0,y2,z2}));
        mol2.getChildList().getAtom(1).getPosition().E(space.makeVector(new double[]{ 0,y2-1.161,z2}));
        mol2.getChildList().getAtom(2).getPosition().E(space.makeVector(new double[]{ 0,y2+1.161,z2 }));
        
//        MoleculeActionTranslateTo translator = new MoleculeActionTranslateTo(space);
//        translator.setDestination(space.makeVector(new double[]{x,0,z}));
//        translator.actionPerformed(mol1);
        PNCO2GCPM p = new PNCO2GCPM(space);
        p.setBox(box);
        p.setComponent(Component.INDUCTION);
        IMoleculeList molecules = box.getMoleculeList();
        P3GCPMAxilrodTeller p3 = p.makeAxilrodTeller();
        double u3 = p3.energy(molecules);
        System.out.println(u3);
        
 
        sim = new Simulation(space);
        SpeciesSpheresRotating species2CO2 = new SpeciesSpheresRotating(space, new ElementSimple("CO2", Carbon.INSTANCE.getMass()+2*Oxygen.INSTANCE.getMass()));
        sim.addSpecies(species2CO2);
        box = new etomica.box.Box(space);
        sim.addBox(box);
        box.setNMolecules(species2CO2, 3);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{100,100,100}));
        IAtomList pair = box.getLeafList();
        IAtomOriented atom0 = (IAtomOriented)pair.getAtom(0);
        IAtomOriented atom1 = (IAtomOriented)pair.getAtom(1);
        IAtomOriented atom2 = (IAtomOriented)pair.getAtom(2);
        atom1.getPosition().E(space.makeVector(new double[]{0,0,z1}));
        atom2.getPosition().E(space.makeVector(new double[]{0,y2,0}));
//        ((IAtomOriented)pair.getAtom(0)).getOrientation().setDirection(space.makeVector(new double[]{Math.cos(22.5/180.0*Math.PI), Math.sin(22.5/180.0*Math.PI),0}));
//        IVectorMutable o1 = space.makeVector(new double[]{-1,0,0});
//        atom1.getOrientation().setDirection(o1);

        AtomTypeAgentManager paramsManagerATM = new AtomTypeAgentManager(null);
        P3AxilrodTeller p3ATM0 = new P3AxilrodTeller(space, paramsManagerATM);
        p3ATM0.setBox(box);
        double alphaCO2 = 2.913;
        paramsManagerATM.setAgent(species2CO2.getLeafType(), new P3AxilrodTeller.MyAgent(alphaCO2, ElectronVolt.UNIT.toSim(13.7)));

        System.out.println(p3ATM0.energy(pair));
        
    }

}
