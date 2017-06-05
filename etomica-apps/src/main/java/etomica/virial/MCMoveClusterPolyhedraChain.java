/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.api.*;
import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomPair;
import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.space.Vector;

public class MCMoveClusterPolyhedraChain extends MCMoveAtom {

    public MCMoveClusterPolyhedraChain(IRandom random, Space _space, double sigma, IPotentialAtomic p2, double[][] uValues) {
        super(random, null, _space);
        this.sigma = sigma;
        dr = space.makeVector();
        this.p2 = p2;
        pair = new AtomPair();
        this.uValues = uValues;
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
        int n = leafAtoms.getAtomCount();
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
        leafAtoms.getAtom(seq[0]).getPosition().E(0);

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                uValues[i][j] = Double.NaN;
            }
        }
        for (int i=1; i<n; i++) {
            pair.atom0 = leafAtoms.getAtom(seq[i-1]);
            pair.atom1 = leafAtoms.getAtom(seq[i]);
            Vector pos = leafAtoms.getAtom(seq[i]).getPosition();
            Vector q = ((AtomOrientedQuaternion)leafAtoms.getAtom(seq[i])).getQuaternion();

            while (true) {
                pos.setRandomInSphere(random);
                pos.TE(sigma);
                pos.PE(leafAtoms.getAtom(seq[i-1]).getPosition());

                randomOrientation(q);
                if (p2.energy(pair) == Double.POSITIVE_INFINITY) break;
            }
            uValues[seq[i-1]][seq[i]] = uValues[seq[i]][seq[i-1]] = Double.POSITIVE_INFINITY;
        }

		((BoxCluster)box).trialNotify();
		return true;
	}
	
    public double getA() {
        return 1;
    }

    public double getB() {
    	return 0.0;
    }
    
    public void rejectNotify() {
        throw new RuntimeException("nope");
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    protected final double sigma;
    protected final Vector dr;
    protected int[] seq;
    protected IPotentialAtomic p2;
    protected final AtomPair pair;
    protected final double[][] uValues;
}
