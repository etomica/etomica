package etomica.modules.glass;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.compute.PotentialCallback;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCallbackAtomStress implements PotentialCallback, AtomStressSource {

    private final Box box;
    private boolean velIncluded = false;
    private final double[][][] stress;
    private final Tensor stressPair;

    public PotentialCallbackAtomStress(Box box) {
        this.box = box;
        int n = box.getLeafList().size();
        int D = box.getSpace().D();
        stress = new double[n][D][D];
        stressPair = box.getSpace().makeTensor();
    }

    public void reset() {
        final int l = stress[0].length;
        for (int i = 0; i < stress.length; i++) {
            for (int j = 0; j < l; j++) {
                for (int k = 0; k < l; k++) {
                    stress[i][j][k] = 0;
                }
            }
        }
    }

    @Override
    public double[][][] getStress() {
        if (!velIncluded) {
            IAtomList atoms = box.getLeafList();
            int D = stress[0].length;
            for (int i = 0; i < atoms.size(); i++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(i);
                Vector velocity = a.getVelocity();
                stressPair.Ev1v2(velocity, velocity);
                stressPair.TE(a.getType().getMass());
                for (int j = 0; j < D; j++) {
                    for (int k = 0; k < D; k++) {
                        stress[i][j][k] = stress[i][j][k] / 2 + stressPair.component(j, k);
                    }
                }
            }
            velIncluded = true;
        }
        return stress;
    }

    @Override
    public void pairCompute(int i, int j, Vector dr, double[] u012) {
        Vector f1 = box.getSpace().makeVector();
        f1.Ea1Tv1(-u012[1] / dr.squared(), dr);
        Tensor t = box.getSpace().makeTensor();
        t.Ev1v2(f1, dr);
        int D = dr.getD();
        for (int k1 = 0; k1 < D; k1++) {
            for (int k2 = 0; k2 < D; k2++) {
                stress[i][k1][k2] += stressPair.component(k1, k2);
                stress[j][k1][k2] += stressPair.component(k1, k2);
            }
        }

    }
}
