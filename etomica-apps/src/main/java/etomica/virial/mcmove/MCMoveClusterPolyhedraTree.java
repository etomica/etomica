/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.mcmove;

import etomica.atom.AtomOrientedQuaternion;
import etomica.atom.AtomPair;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.Potential2Soft;
import etomica.space.Vector;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

public class MCMoveClusterPolyhedraTree extends MCMoveBox {

    protected final IRandom random;
    protected final double sigma;
    protected int[][] bonds;
    protected int[] degree, a;
    protected int[] inserted;
    protected Potential2Soft p2;
    protected final AtomPair pair;
    protected final double[][] uValues;

    public MCMoveClusterPolyhedraTree(IRandom random, Box box, double sigma, Potential2Soft p2, double[][] uValues) {
        super();
        this.random = random;
        this.sigma = sigma;
        this.p2 = p2;
        pair = new AtomPair();
        this.uValues = uValues;
        setBox(box);
    }
    
    public void setBox(Box box) {
        super.setBox(box);
        int n = box.getLeafList().size();
        degree = new int[n];
        a = new int[n-2];
        inserted = new int[n];
        bonds = new int[n*(n-1)/2][2];
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
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

        leafAtoms.get(0).getPosition().E(0);
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
                IAtom atom1 = leafAtoms.get(nbr);
                IAtom atom2 = leafAtoms.get(nbr2);
                // insert nbr2 around nbr
                Vector q = ((AtomOrientedQuaternion)leafAtoms.get(nbr2)).getQuaternion();
                Vector pos = leafAtoms.get(nbr2).getPosition();

                Vector dr = box.getSpace().makeVector();
                while (true) {
                    dr.setRandomInSphere(random);
                    dr.TE(sigma);
                    pos.Ev1Pv2(atom1.getPosition(), dr);

                    randomOrientation(q);
                    if (p2.u(dr, atom1, atom2) == Double.POSITIVE_INFINITY) break;
                }
                uValues[nbr2][nbr] = uValues[nbr][nbr2] = Double.POSITIVE_INFINITY;
                
                inserted[numInserted] = nbr2;
                numInserted++;
            }
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
}
