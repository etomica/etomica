/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomPair;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.IPotentialAtomic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveClusterPolyhedraTree extends MCMoveAtom {

    public MCMoveClusterPolyhedraTree(IRandom random, Space _space, double sigma, IPotentialAtomic p2, double[][] uValues) {
        super(random, null, _space);
        this.sigma = sigma;
        this.p2 = p2;
        pair = new AtomPair();
        this.uValues = uValues;
    }
    
    public void setBox(Box box) {
        super.setBox(box);
        int n = box.getLeafList().getAtomCount();
        degree = new int[n];
        a = new int[n-2];
        inserted = new int[n];
        bonds = new int[n*(n-1)/2][2];
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
        for (int i=0; i<n; i++) {
            degree[i] = 1;
        }
        for (int i=0; i<n-2; i++) {
            a[i] = random.nextInt(n);
            degree[a[i]]++;
        }

        int numBonds = 0;
        for (int i=0; i<n-2; i++) {
            int ii = a[i];
            for (int j=0; j<n; j++) {
                if (degree[j] == 1) {
                    bonds[numBonds][0] = ii;
                    bonds[numBonds][1] = j;
                    numBonds++;
                    degree[ii]--;
                    degree[j]--;
                    break;
                }
            }
        }
        int u = -1, v = -1;
        for (int i=0; i<n; i++) {
            if (degree[i] == 1) {
                if (u==-1) {
                    u = i;
                }
                else {
                    v = i;
                }
            }
        }
        bonds[numBonds][0] = u;
        bonds[numBonds][1] = v;
        numBonds++;

        leafAtoms.getAtom(0).getPosition().E(0);
        inserted[0] = 0;
        int numInserted = 1;
        // inserted is a list of points that have inserted, but not coordinated
        // coordinated is a list of points that have been inserted and coordinated
        int coordinatedMask = 0;
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                uValues[i][j] = Double.NaN;
            }
        }
        while (numInserted > 0) {
            int nbr = inserted[numInserted-1];
            numInserted--;
            coordinatedMask |= 1<<nbr;
            for (int iBond = 0; iBond<numBonds; iBond++) {
                int[] b = bonds[iBond];
                int nbr2 = -1;
                if (b[0] == nbr) {
                    nbr2 = b[1];
                }
                else if (b[1] == nbr) {
                    nbr2 = b[0];
                }
                else {
                    continue;
                }
                if ((coordinatedMask & (1<<nbr2)) != 0) {
                    // already inserted nbr2, move along
                    continue;
                }
                pair.atom0 = leafAtoms.getAtom(nbr);
                pair.atom1 = leafAtoms.getAtom(nbr2);
                // insert nbr2 around nbr
                Vector q = ((AtomOrientedQuaternion)leafAtoms.getAtom(nbr2)).getQuaternion();
                Vector pos = leafAtoms.getAtom(nbr2).getPosition();

                while (true) {
                    pos.setRandomInSphere(random);
                    pos.TE(sigma);
                    pos.PE(leafAtoms.getAtom(nbr).getPosition());

                    randomOrientation(q);
                    double energy = p2.energy(pair);
                    if (energy == Double.POSITIVE_INFINITY) break;
                }
                uValues[nbr2][nbr] = uValues[nbr][nbr2] = Double.POSITIVE_INFINITY;
                
                inserted[numInserted] = nbr2;
                numInserted++;
            }
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
    protected int[][] bonds;
    protected int[] degree, a;
    protected int[] inserted;
    protected IPotentialAtomic p2;
    protected final AtomPair pair;
    protected final double[][] uValues;
}
