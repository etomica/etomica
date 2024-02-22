/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveHO extends MCMoveBox {
    protected int nBeads;
    protected double temperature, omega2;
    protected final IRandom random;
    double[] lambda;
    double[][] eigenvectors, eigenvectorsInv;
    PotentialCompute pm;
    protected final Vector[] oldPositions;
    protected double uOld;
    protected double uaOld = Double.NaN;
    protected double uaNew = Double.NaN;
    protected double duTotal;
    protected double mass, beta, omegaN;
    protected Box box;


    public MCMoveHO(Space space, PotentialCompute pm, IRandom random, double temperature, double omega2, Box box, double hbar) {
        super();
        this.pm = pm;
        this.random = random;
        this.temperature = temperature;
        this.omega2 = omega2;
        this.box = box;
        nBeads = this.box.getMoleculeList().get(0).getChildList().size();
        oldPositions = new Vector[nBeads];
        for (int i=0; i < nBeads; i++){
            oldPositions[i] = space.makeVector();
        }

        lambda = new double[nBeads];
        beta = 1.0/temperature;
        omegaN = Math.sqrt(nBeads)/(hbar*beta);
        mass = box.getLeafList().get(0).getType().getMass();// m/n
        int nK = nBeads/2;
        lambda[nK] = mass*omega2; //k=0
        double lambda_k;
        for(int k = 1; k <= (nBeads-1)/2; k++){
            lambda_k = 4.0*mass*nBeads*omegaN*omegaN*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda[nK-k] = lambda_k;
            lambda[nK+k] = lambda_k;
        }
        if (nBeads % 2 == 0){
            int k = nK;
            lambda_k = 4.0*mass*nBeads*omegaN*omegaN*Math.sin(Math.PI*k/nBeads)*Math.sin(Math.PI*k/nBeads) + mass*omega2;
            lambda[0] = lambda_k;
        }

        eigenvectors = new double[nBeads][nBeads];
        eigenvectorsInv = new double[nBeads][nBeads];//inv=transpose
        for (int i = 0; i < nBeads; i++) {
            eigenvectors[i][nK] = 1.0/Math.sqrt(nBeads);//k=0
            eigenvectorsInv[nK][i] = 1.0/Math.sqrt(nBeads);//k=0
            for (int k = 1; k <= (nBeads-1)/2; k++) {
                double arg = 2.0*Math.PI/nBeads*i*k;
                eigenvectors[i][nK-k] = Math.sqrt(2.0)*Math.sin(-arg)/Math.sqrt(nBeads);
                eigenvectors[i][nK+k] = Math.sqrt(2.0)*Math.cos(arg)/Math.sqrt(nBeads);
                eigenvectorsInv[nK-k][i] = Math.sqrt(2.0)*Math.sin(-arg)/Math.sqrt(nBeads); //sin then cos .. it's a must!
                eigenvectorsInv[nK+k][i] = Math.sqrt(2.0)*Math.cos(arg)/Math.sqrt(nBeads);
            }
            if (nBeads % 2 == 0){ //even
                int k = nK;
                double arg = 2.0*Math.PI/nBeads*i*k;
                eigenvectors[i][0] =  Math.pow(-1, i)/Math.sqrt(nBeads);// Math.pow(-1, i)=Math.cos(arg);
                eigenvectorsInv[0][i] =  Math.pow(-1, i)/Math.sqrt(nBeads);
            }
        }
    }

    public boolean doTrial() {
        double uOld = pm.getLastEnergy();
        IAtomList atomList = box.getLeafList();
        IAtom atomi, atomip1;
        Vector ri, rip1;
        double uhOld = 0;
        for (int i = 0; i < nBeads; i++) {
            atomi = atomList.get(i);
            atomip1 = (i == nBeads - 1) ? atomList.get(0) : atomList.get(i + 1);
            ri = atomi.getPosition();
            oldPositions[i].E(ri);
            rip1 = atomip1.getPosition();
            Vector dr = box.getSpace().makeVector();
            dr.Ev1Mv2(ri, rip1);
            uhOld += 0.5*mass*nBeads*(omegaN*omegaN*dr.squared() + omega2*ri.squared()/nBeads);
        }
        uaOld = uOld - uhOld;

        double uhNew = 0;
        Vector[] q = box.getSpace().makeVectorArray(nBeads);

        int nK = nBeads/2;
        double fack = 1.0/Math.sqrt(beta*lambda[nK]);

        for (int l=0; l<q[0].getD(); l++) {
            q[nK].setX(l, fack * random.nextGaussian());
        }

        if (omega2 == 0){
            q[nK].E(0);
            for (int i = 0; i < nBeads; i++) {
                atomi = atomList.get(i);
                ri = atomi.getPosition();
                q[nK].PEa1Tv1(1/Math.sqrt(nBeads), ri);
            }
        }

        uhNew += 1.0/2.0*lambda[nK]*q[nK].dot(q[nK]);
        for (int k = 1; k <= (nBeads-1)/2; k++) {
            fack = 1.0/Math.sqrt(beta*lambda[nK-k]);
            for (int l=0; l<q[0].getD(); l++) {
                q[nK-k].setX(l, fack * random.nextGaussian());
            }
            uhNew += 1.0/2.0*lambda[nK-k]*q[nK-k].dot(q[nK-k]);
            fack = 1.0/Math.sqrt(beta*lambda[nK+k]);
            for (int l=0; l<q[0].getD(); l++) {
                q[nK+k].setX(l, fack * random.nextGaussian());
            }
            uhNew += 1.0/2.0*lambda[nK+k]*q[nK+k].dot(q[nK+k]);
        }
        if (nBeads % 2 == 0){ //even
            fack = 1.0/Math.sqrt(beta*lambda[0]);
            for (int l=0; l<q[0].getD(); l++) {
                q[0].setX(l, fack * random.nextGaussian());
            }
            uhNew += 1.0/2.0*lambda[0]*q[0].dot(q[0]);
        }

        for (int i = 0; i < nBeads; i++) {
            atomi = atomList.get(i);
            ri = atomi.getPosition();
            ri.E(0);
            for (int k = 0; k < nBeads; k++) {
                ri.PEa1Tv1(eigenvectors[i][k], q[k]);
            }
        }

        boolean debug = false;
        if (debug){
            double uReal = 0;
            for (int i = 0; i < nBeads; i++) {
                atomi = atomList.get(i);
                atomip1 = (i == nBeads - 1) ? atomList.get(0) : atomList.get(i + 1);
                ri = atomi.getPosition();
                oldPositions[i].E(ri);
                rip1 = atomip1.getPosition();
                Vector dr = box.getSpace().makeVector();
                dr.Ev1Mv2(ri, rip1);
                uReal += 0.5*mass*nBeads*(omegaN*omegaN*dr.squared() + omega2*ri.squared()/nBeads);
            }
            double uNM = uhNew;
            double foo = uNM/uReal;
            System.out.println(foo );
            System.out.println();
        }

        for (int k = 0; k < nBeads; k++) {
            pm.updateAtom(atomList.get(k));
        }

        double uNew = pm.computeAll(false);
        uaNew = uNew - uhNew;

        duTotal = uNew - uOld;

        return true;
    }

    public double energyChange() {return duTotal;}

    public double getChi(double temperature) {
        return Math.exp(-(uaNew - uaOld) / temperature);
    }

    public void acceptNotify() { /* do nothing */}

    public void rejectNotify() {
        IAtomList atomList = box.getLeafList();
        IAtom atomi;
        for (int i = 0; i < nBeads; i++) {
            atomi = atomList.get(i);
            atomi.getPosition().E(oldPositions[i]);
        }
        pm.init();
        pm.computeAll(false);
    }

}
