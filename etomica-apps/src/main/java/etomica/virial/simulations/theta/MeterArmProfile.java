
package etomica.virial.simulations.theta;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.space.Vector;
import etomica.units.dimensions.Length;

/**
 * Meter that computes the average radial distance of each monomer
 * along an arm from the core center (atom 0).
 *
 * For a star polymer with numArms arms each of length armLength:
 *   - atom 0                        = core
 *   - atom 1 + i*armLength + (k-1)  = monomer k (1-indexed) of arm i (0-indexed)
 *
 * getData() returns an array of length armLength.
 * Entry k-1 (0-indexed) = average distance from core of the k-th monomer,
 * averaged over all numArms arms in the current configuration.
 *
 * Wire this like MeterRadiusGyration:
 *   MeterArmProfile meterProfile = new MeterArmProfile(sim.box(), numArms, armLength);
 *   AccumulatorAverageFixed accProfile = new AccumulatorAverageFixed(steps/1000);
 *   DataPumpListener pumpProfile = new DataPumpListener(meterProfile, accProfile, 10);
 *   integrator.getEventManager().addListener(pumpProfile);
 */
public class MeterArmProfile implements IDataSource {

    private final Box box;
    private final int numArms;
    private final int armLength;

    // data array of length armLength: index k-1 holds avg r for monomer k
    private final DataDoubleArray data;
    private final DataDoubleArray.DataInfoDoubleArray dataInfo;
    private final DataTag tag;

    public MeterArmProfile(Box box, int numArms, int armLength) {
        this.box = box;
        this.numArms = numArms;
        this.armLength = armLength;

        // one slot per bead position along the arm (1 .. armLength)
        data = new DataDoubleArray(armLength);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray(
                "ArmProfile", Length.DIMENSION, new int[]{armLength});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    /**
     * Called every interval (e.g. every 10 steps via DataPumpListener).
     * Returns array of length armLength: entry [k-1] = average distance
     * from core of the k-th monomer along the arm, averaged over all arms.
     */
    @Override
    public IData getData() {
        IAtomList atoms = box.getLeafList();
        // core is always atom 0
        Vector corePos = atoms.get(0).getPosition();

        double[] y = data.getData();
        // zero out before accumulating
        for (int k = 0; k < armLength; k++) {
            y[k] = 0.0;
        }

        // sum distances over all arms
        for (int i = 0; i < numArms; i++) {
            for (int k = 1; k <= armLength; k++) {
                // leaf atom index for arm i, monomer k (1-indexed from core)
                int atomIndex = 1 + i * armLength + (k - 1);
                Vector pos = atoms.get(atomIndex).getPosition();
                // distance from core center
                double r = Math.sqrt(pos.Mv1Squared(corePos));
                y[k - 1] += r;
            }
        }

        // average over all arms
        for (int k = 0; k < armLength; k++) {
            y[k] /= numArms;
        }

        return data;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }
}