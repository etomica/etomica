package etomica.spin.heisenberg_interacting.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorRigidIterative;
import etomica.integrator.IntegratorRigidIterative.MoleculeAgent;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.units.dimensions.Null;
import sun.security.krb5.KdcComm;

public class MeterMappedAveragingCorrelation implements IDataSource, AtomLeafAgentManager.AgentSource<IntegratorRigidIterative.MoleculeAgent> {
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected double bJ;
    protected double bt;
    protected final int L;
    protected final int N;
    protected double temperature;
    private Box box;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected PotentialCalculationTorqueSum torqueSum;
    protected final IteratorDirective allAtoms;
    protected final int[] offset;
    protected final boolean formula1;
    protected AtomLeafAgentManager<IntegratorRigidIterative.MoleculeAgent> leafAgentManager;


    public MeterMappedAveragingCorrelation(Simulation sim, double temperature, double interactionS, PotentialMaster potentialMaster, boolean formula1) {
        this.box = sim.getBox(0);
        this.space = sim.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        this.formula1 = formula1;
        bt = 1 / temperature;
        bJ = interactionS * bt;
        N = box.getLeafList().getAtomCount();
        L = (int) Math.round(Math.sqrt(N));
        int distance = L / 2 + 1;
        data = new DataDoubleArray(-1 + (distance + 1) * distance / 2);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("CIJ", Null.DIMENSION, new int[]{-1 + (distance + 1) * distance / 2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        allAtoms = new IteratorDirective();
        offset = new int[distance];
        for (int i = 0, o = -1; i < distance; i++) {
            offset[i] = o - i;
            o += (distance - i);
        }

        leafAgentManager = new AtomLeafAgentManager<IntegratorRigidIterative.MoleculeAgent>(this, box, IntegratorRigidIterative.MoleculeAgent.class);
        torqueSum = new PotentialCalculationTorqueSum();
        torqueSum.setAgentManager(leafAgentManager);
    }

    public IData getData() {
        double[] x = data.getData();
        int distance = L / 2 + 1;
        torqueSum.reset();
        potentialMaster.calculate(box, allAtoms, torqueSum);
        IAtomList leafList = box.getLeafList();
        for (int i = 0; i < x.length; i++) {
            x[i] = 0;
        }
        for (int j = 0; j < N; j++) {
            int jRow = j / L, jCol = j % L;
            IntegratorRigidIterative.MoleculeAgent agentAtomJ = leafAgentManager.getAgent(leafList.getAtom(j));
            double fj = bt * agentAtomJ.torque().getX(0);
            double costj = ((IAtomOriented) leafList.getAtom(j)).getOrientation().getDirection().getX(0);
            double sintj = ((IAtomOriented) leafList.getAtom(j)).getOrientation().getDirection().getX(1);

            for (int k = 0; k < j; k++) {
                int kRow = k / L, kCol = k % L;
                IntegratorRigidIterative.MoleculeAgent agentAtomK = leafAgentManager.getAgent(leafList.getAtom(k));
                double fk = bt * agentAtomK.torque().getX(0);
                double costk = ((IAtomOriented) leafList.getAtom(k)).getOrientation().getDirection().getX(0);
                double sintk = ((IAtomOriented) leafList.getAtom(k)).getOrientation().getDirection().getX(1);
                int JMKRow = jRow - kRow;
                int JMKCol = jCol - kCol;
                if (JMKRow > distance) JMKRow = L - JMKRow;
                if (JMKCol > distance) {
                    JMKCol = L - JMKCol;
                } else if (JMKCol < -distance) {
                    JMKCol += L;
                } else if (JMKCol < 0) {
                    JMKCol *= -1;
                }

                if (JMKCol < JMKRow) {
                    int tmp = JMKCol;
                    JMKCol = JMKRow;
                    JMKRow = tmp;
                }
                int index = offset[JMKRow] + JMKCol;
                if (formula1) {
                    x[index] += 0.5 * (fj - fk) * (sintj * costk - costj * sintk);
                } else {
                    x[index] -= fj * fk * (sintj * sintk + costj * costk);
                }

            }
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public IntegratorRigidIterative.MoleculeAgent makeAgent(IAtom a, Box agentBox) {

        return new IntegratorRigidIterative.MoleculeAgent(space);
    }

    public void releaseAgent(IntegratorRigidIterative.MoleculeAgent agent, IAtom atom, Box agentBox) {

    }
}
