/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.Potential2Soft;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

public class MCMoveClusterPolyhedraChain extends MCMoveBox {

    protected final IRandom random;
    protected final double sigma;
    protected final Vector dr;
    protected int[] seq;
    protected Potential2Soft p2;
    protected final double[][] uValues;

    public MCMoveClusterPolyhedraChain(IRandom random, Box box, double sigma, Potential2Soft p2, double[][] uValues) {
        super();
        this.random = random;
        this.sigma = sigma;
        dr = box.getSpace().makeVector();
        this.p2 = p2;
        this.uValues = uValues;
        setBox(box);
    }

    protected void randomOrientation(Vector q) {
        double u1 = random.nextDouble();
        double u2 = 2*Math.PI*random.nextDouble();
        double u3 = 2*Math.PI*random.nextDouble();
        double s1 = Math.sqrt(u1);
        double s2 = Math.sqrt(1-u1);
        q.setX(0, s1*Math.sin(u2));
        q.setX(1, s1*Math.cos(u2));
        q.setX(2, s2*Math.sin(u3));
        q.setX(3, s2*Math.cos(u3));
    }
    
    public boolean doTrial() {
        
        IAtomList leafAtoms = box.getLeafList();
        int n = leafAtoms.size();
        if (seq == null) {
            seq = new int[n];
        }
        for (int i=0; i<n; i++) {
            seq[i] = i;
        }
        for (int i=0; i<n; i++) {
            int j = i+random.nextInt(n-i);
            int k = seq[j];
            seq[j] = seq[i];
            seq[i] = k;
        }
        leafAtoms.get(seq[0]).getPosition().E(0);

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                uValues[i][j] = Double.NaN;
            }
        }
        Vector dr = box.getSpace().makeVector();
        for (int i=1; i<n; i++) {
            IAtom atom1 = leafAtoms.get(seq[i-1]);
            IAtom atom2 = leafAtoms.get(seq[i]);
            Vector pos = leafAtoms.get(seq[i]).getPosition();
            Vector q = ((AtomOrientedQuaternion)leafAtoms.get(seq[i])).getQuaternion();

            while (true) {
                dr.setRandomInSphere(random);
                dr.TE(sigma);
                pos.Ev1Pv2(atom1.getPosition(), dr);

                randomOrientation(q);
                if (p2.u(dr, atom1, atom2) == Double.POSITIVE_INFINITY) break;
            }
            uValues[seq[i-1]][seq[i]] = uValues[seq[i]][seq[i-1]] = Double.POSITIVE_INFINITY;
        }

		((BoxCluster)box).trialNotify();
		return true;
	}

    public double getChi(double temperature) {
        return 1;
    }

    public void rejectNotify() {
        throw new RuntimeException("nope");
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }
}
