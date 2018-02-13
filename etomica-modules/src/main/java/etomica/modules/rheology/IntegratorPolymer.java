/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.rheology;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.IntegratorMD;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 * Integrator for Brownian dynamics of a polymer in a flow field.
 *
 * @author Andrew Schultz
 */
public class IntegratorPolymer extends IntegratorMD {

    public IntegratorPolymer(PotentialMaster potentialMaster, IRandom random,
                             double timeStep, double temperature, Space _space, Box box) {
        super(potentialMaster, random, timeStep, temperature, _space, box);
        center = _space.makeVector();
        drPrev = _space.makeVector();
        dr = _space.makeVector();
        ds = _space.makeVector();
        s = new Vector[0];
        r = new Vector[0];
        W = new Vector[0];
        fQ = new double[0];
    }

    /* no point in trying to thermostat */
    protected void doThermostatInternal() {}

    protected void doStepInternal() {
        super.doStepInternal();

        double z = 1;
        if (a > 0) {
            z = (1-a*a)/(1+a);
        }
        double sqa = a>0 ? Math.sqrt(a) : 0;
        double srdt = shearRate*timeStep;
        if (sqa > 1.0/160.0) {
            srdt /= 162*sqa;
        }
        double srdt2 = 0.5*srdt;
        IMoleculeList molecules = box.getMoleculeList();
        for (int i=0; i<molecules.getMoleculeCount(); i++) {
            center.E(0);
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atoms = molecule.getChildList();
            if (s.length != atoms.getAtomCount()) {
                s = new Vector[atoms.getAtomCount()];
                r = new Vector[atoms.getAtomCount()];
                W = new Vector[atoms.getAtomCount()];
                for (int j=0; j<atoms.getAtomCount(); j++) {
                    s[j] = space.makeVector();
                    r[j] = space.makeVector();
                    W[j] = space.makeVector();
                }
                fQ = new double[atoms.getAtomCount()];
            }
            for (int j=0; j<atoms.getAtomCount(); j++) {
                s[j].E(atoms.getAtom(j).getPosition());
                for (int k=0; k<3; k++) {
                    W[j].setX(k, sqdt*random.nextGaussian());
                }
            }
            //calculate 0-1 bond vector before moving 0 and 1.

            // predictor step
            for (int j=0; j<atoms.getAtomCount(); j++) {

                if (a < 0) {
                    r[j].setX(0, srdt*s[j].getX(1));
                    r[j].setX(1, a*srdt*s[j].getX(0));
                }
                else {
                    r[j].setX(0, srdt*(z*s[j].getX(1) + sqa*s[j].getX(0)));
                    r[j].setX(1, -sqa*srdt*s[j].getX(1));
                }
                r[j].setX(2,0);

                r[j].PE(s[j]);

                if (j > 0) {
                    r[j].PEa1Tv1(fQ[j-1], ds);
                }

                if (j+1 < atoms.getAtomCount()) {
                    ds.Ev1Mv2(s[j+1], s[j]);
                    fQ[j] = (1 + b*ds.squared())*omdth;
                    r[j].PEa1Tv1(-fQ[j], ds);
                }
                
                r[j].PE(W[j]);
            }
            
            // corrector step
            double fR = 0;
            for (int j=0; j<atoms.getAtomCount(); j++) {
                Vector q = atoms.getAtom(j).getPosition();
                if (a < 0) {
                    q.setX(0, srdt2*(s[j].getX(1)+r[j].getX(1)));
                    q.setX(1, a*srdt2*(s[j].getX(0)+r[j].getX(0)));
                }
                else {
                    q.setX(0, srdt2*(z*(s[j].getX(1)+r[j].getX(1)) + sqa*(s[j].getX(0)+r[j].getX(0))));
                    q.setX(1, -sqa*srdt2*(s[j].getX(1)+r[j].getX(1)));
                }
                q.setX(2, 0);
                q.PE(s[j]);
                if (j > 0) {
                    q.PEa1Tv1(fQ[j-1], ds);
                    q.PEa1Tv1(fR, dr);
                }
                if (j+1 < atoms.getAtomCount()) {
                    ds.Ev1Mv2(s[j+1], s[j]);
                    fQ[j] *= 0.5;
                    q.PEa1Tv1(-fQ[j], ds);
                    dr.Ev1Mv2(s[j+1], s[j]);
                    fR = (1 + b*dr.squared())*omdth*0.5;
                    q.PEa1Tv1(-fR, dr);
                }
                q.PE(W[j]);
                center.PE(q);
            }

            // maintain center at 0
            center.TE(-1.0/atoms.getAtomCount());

            for (int j=0; j<atoms.getAtomCount(); j++) {
                atoms.getAtom(j).getPosition().PE(center);
            }
        }
        
    }

    public void setTimeStep(double newTimeStep) {
        super.setTimeStep(newTimeStep);
        omdth = -0.25*timeStep;
        sqdt = Math.sqrt(0.5*timeStep);
    }

    public void setShearRateNumber(double newShearRate) {
        shearRate = newShearRate;
    }
    
    public double getShearRate() {
        if (a<0) return shearRate;
        double sqa = Math.sqrt(a);
        double sr = shearRate;
        if (sqa > 1.0/160.0) {
            sr /= 162*sqa;
        }
        return sr;
    }

    public double getShearRateNumber() {
        return shearRate;
    }

    public double getA() {
        return a;
    }

    public void setA(double newA) {
        a = newA;
    }

    public double getB() {
        return b;
    }

    public void setB(double newB) {
        b = newB;
    }

    private static final long serialVersionUID = 1L;
    protected final Vector drPrev, dr, ds, center;
    protected double omdth, sqdt;
    protected double shearRate;
    protected double a, b;
    protected Vector[] W;
    protected Vector[] s, r;
    protected double[] fQ;
}
