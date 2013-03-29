package etomica.virial;

import java.util.ArrayList;
import java.util.List;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IRandom;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;
import etomica.space.IVectorRandom;

public class MCMoveClusterAtomHSTree extends MCMoveAtom {

    public MCMoveClusterAtomHSTree(IRandom random, ISpace _space, double sigma) {
        super(random, null, _space);
        this.sigma = sigma;
        dr = (IVectorRandom)space.makeVector();
        bonds = new ArrayList<int[]>();
    }
    
    public void setBox(IBox box) {
        super.setBox(box);
        int n = box.getLeafList().getAtomCount();
        degree = new int[n];
        a = new int[n-2];
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

        bonds.clear();
        for (int i=0; i<n-2; i++) {
            int ii = a[i];
            for (int j=0; j<n; j++) {
                if (degree[j] == 1) {
                    bonds.add(new int[]{ii,j});
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
        bonds.add(new int[]{u,v});

        leafAtoms.getAtom(0).getPosition().E(0);
        List<Integer> inserted = new ArrayList<Integer>();
        inserted.add(0);
        // inserted is a list of points that have inserted, but not coordinated
        List<Integer> coordinated = new ArrayList<Integer>();
        // coordinted is a list of points that have been inserted and coordinated
//        int[] c = new int[n];
        while (inserted.size() > 0) {
            int nbr = inserted.remove(0);
            coordinated.add(nbr);
            for (int[] b : bonds) {
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
                if (coordinated.contains(nbr2)) {
                    // already inserted nbr2, move along
                    continue;
                }
                // insert nbr2 around nbr
                IVectorRandom pos = (IVectorRandom)leafAtoms.getAtom(nbr2).getPosition();

                pos.setRandomInSphere(random);
                pos.TE(sigma);
                pos.PE(leafAtoms.getAtom(nbr).getPosition());
                inserted.add(nbr2);
                
//                c[nbr]++;
//                c[nbr2]++;
            }
        }
        
        /*t++;
        for (int i=0; i<n; i++) {
            if (c[i]==n-1) {
                System.out.println("b "+i);
                bb[i]++;
                br++;
            }
        } */
//        if (t%1000 == 0) {
//            System.out.print(String.format("%d  %5.3f %5.3f %5.3f %5.3f %5.3f\n", t, t/((double)br), t/((double)bb[0]), t/((double)bb[1]), t/((double)bb[2]), t/((double)bb[3])));
//        }

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
    protected final IVectorRandom dr;
    protected final List<int[]> bonds;
    protected int[] degree, a;
    /*public static long br, t;
    public static long[] bb = new long[5]; */
}
