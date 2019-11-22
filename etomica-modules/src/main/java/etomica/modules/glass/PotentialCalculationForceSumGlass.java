package etomica.modules.glass;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculationForcePressureSum;
import etomica.potential.PotentialSoft;
import etomica.space.Tensor;
import etomica.space.Vector;

public class PotentialCalculationForceSumGlass extends PotentialCalculationForcePressureSum implements AtomStressSource {

    private final Box box;
    private boolean velIncluded = false;
    private final double[][] stress;
    private final Tensor stressPair;

    public PotentialCalculationForceSumGlass(Box box) {
        super(box.getSpace());
        this.box = box;
        int n = box.getLeafList().size();
        int m = box.getSpace().D() == 3 ? 3 : 1;
        stress = new double[n][m];
        stressPair = box.getSpace().makeTensor();
    }

    public void reset() {
        super.reset();
        final int l = stress[0].length;
        for (int i = 0; i < stress.length; i++) {
            for (int j = 0; j < l; j++) {
                stress[i][j] = 0;
            }
        }
    }

    /**
     * Adds forces due to given potential acting on the atoms produced by the iterator.
     * Implemented for only 1- and 2-body potentials.
     */
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        if (!(potential instanceof PotentialSoft)) return;
        velIncluded = false;
        PotentialSoft potentialSoft = (PotentialSoft) potential;
        stressPair.E(0);
        Vector[] f = potentialSoft.gradient(atoms, stressPair);
        pressureTensor.PE(stressPair);
        integratorAgentManager.getAgent(atoms.get(0)).ME(f[0]);
        integratorAgentManager.getAgent(atoms.get(1)).ME(f[1]);
        int idx0 = atoms.get(0).getLeafIndex();
        int idx1 = atoms.get(1).getLeafIndex();
        stress[idx0][0] += stressPair.component(0, 1);
        stress[idx1][0] += stressPair.component(0, 1);
        if (f.length == 3) {
            stress[idx0][1] += stressPair.component(0, 2);
            stress[idx1][1] += stressPair.component(0, 2);
            stress[idx0][2] += stressPair.component(1, 2);
            stress[idx1][2] += stressPair.component(1, 2);
        }
    }

    @Override
    public double[][] getStress() {
        if (!velIncluded) {
            IAtomList atoms = box.getLeafList();
            for (int i = 0; i < atoms.size(); i++) {
                IAtomKinetic a = (IAtomKinetic) atoms.get(i);
                Vector velocity = a.getVelocity();
                stressPair.Ev1v2(velocity, velocity);
                stressPair.TE(a.getType().getMass());
                stress[i][0] += stressPair.component(0, 1);
                if (stress[i].length == 3) {
                    stress[i][1] += stressPair.component(0, 2);
                    stress[i][2] += stressPair.component(1, 2);
                }
            }
            velIncluded = true;
        }
        return stress;
    }
}
