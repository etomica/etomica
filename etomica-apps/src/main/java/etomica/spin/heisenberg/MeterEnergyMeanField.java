package etomica.spin.heisenberg;

import etomica.atom.IAtomOriented;
import etomica.space.Vector;

//public class MeterEnergyMeanField implements IDataSource, AtomLeafAgentManager.AgentSource<MeterEnergyMeanField.ForceTorque> {
public class MeterEnergyMeanField {

    /**
     * Used to compute and store sums of cos(theta) and sin(theta) over neighbors of each atom.
     * Needed for mapped-averaging calculation of heat capacity.
     */
    public static class PotentialCalculationCSsum {
        protected final double[] nbrCsum, nbrSsum;
        public PotentialCalculationCSsum(double[] nbrCsum, double[] nbrSsum) {
            this.nbrCsum = nbrCsum;
            this.nbrSsum = nbrSsum;
        }
        public void go(IAtomOriented iatom, IAtomOriented jatom) {
            int i = iatom.getLeafIndex();
            int j = jatom.getLeafIndex();
            Vector io = iatom.getOrientation().getDirection();
            Vector jo = jatom.getOrientation().getDirection();
            nbrCsum[i] += jo.getX(0);
            nbrCsum[j] += io.getX(0);
            nbrSsum[i] += jo.getX(1);
            nbrSsum[j] += io.getX(1);
        }

        public void reset() {
            for (int i = 0; i < nbrCsum.length; i++) {
                nbrCsum[i] = 0;
                nbrSsum[i] = 0;
            }
        }
    }
}