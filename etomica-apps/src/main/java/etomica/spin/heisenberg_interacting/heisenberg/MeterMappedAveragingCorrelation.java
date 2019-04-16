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
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space1d.Vector1D;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingCorrelation implements IDataSource, AtomLeafAgentManager.AgentSource<MeterMappedAveragingCorrelation.CorrelationAgent> {
    protected final DataDoubleArray data;
    protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected double bJ;
    protected double bt;
    protected final int L;
    protected final int N;
    protected final int arraySize;
    protected double temperature;
    private Box box;
    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected PotentialCalculationTorqueSum torqueSum;
    protected final IteratorDirective allAtoms;
    protected final int[] offset;
    protected final int formula;
    protected final int[] nPairs;
    protected AtomLeafAgentManager<CorrelationAgent> leafAgentManager;


    public MeterMappedAveragingCorrelation(Simulation sim, double temperature, double interactionS, PotentialMaster potentialMaster, int formula) {
        this.box = sim.getBox(0);
        this.space = sim.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        this.formula = formula;
        bt = 1 / temperature;
        bJ = interactionS * bt;
        N = box.getLeafList().getAtomCount();
        L = (int) Math.round(Math.sqrt(N));
        int distance = L / 2 + 1;
//        nPairs = new int[-1 + (distance + 1) * distance / 2];
//        data = new DataDoubleArray(-1 + (distance + 1) * distance / 2);
//        dataInfo = new DataDoubleArray.DataInfoDoubleArray("CIJ", Null.DIMENSION, new int[]{-1 + (distance + 1) * distance / 2});
        arraySize = -1 + (distance + 1) * distance / 2;
        nPairs = new int[3 * arraySize];
        data = new DataDoubleArray(3 * arraySize);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("CIJ", Null.DIMENSION, new int[]{3 * arraySize});
        tag = new DataTag();
        dataInfo.addTag(tag);
        allAtoms = new IteratorDirective();
        offset = new int[distance];
        for (int i = 0, o = -1; i < distance; i++) {
            offset[i] = o - i;
            o += (distance - i);
        }

        leafAgentManager = new AtomLeafAgentManager<CorrelationAgent>(this, box, CorrelationAgent.class);
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
            nPairs[i] = 0;
        }
        for (int j = 0; j < N; j++) {
            int jRow = j / L, jCol = j % L;
            CorrelationAgent agentAtomJ = leafAgentManager.getAgent(leafList.getAtom(j));
            double fj = bt * agentAtomJ.torque().getX(0);
            double costj = ((IAtomOriented) leafList.getAtom(j)).getOrientation().getDirection().getX(0);
            double sintj = ((IAtomOriented) leafList.getAtom(j)).getOrientation().getDirection().getX(1);

            for (int k = 0; k < j; k++) {
                int kRow = k / L, kCol = k % L;
                CorrelationAgent agentAtomK = leafAgentManager.getAgent(leafList.getAtom(k));
                double fk = bt * agentAtomK.torque().getX(0);
                double costk = ((IAtomOriented) leafList.getAtom(k)).getOrientation().getDirection().getX(0);
                double sintk = ((IAtomOriented) leafList.getAtom(k)).getOrientation().getDirection().getX(1);
                int JMKRow = jRow - kRow;
                int JMKCol = jCol - kCol;
                if (JMKRow >= distance) JMKRow = L - JMKRow;
                if (JMKCol >= distance) {
                    JMKCol = L - JMKCol;
                } else if (JMKCol <= -distance) {
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
                if (formula == 0) {
                    x[index] += costj * costk + sintj * sintk;
                } else if (formula == 1) {
                    x[index] += 0.5 * (fj - fk) * (sintj * costk - costj * sintk);
                } else {
                    x[index] -= fj * fk * (sintj * sintk + costj * costk);
                }
                nPairs[index] += 1;
                x[index + arraySize] = JMKRow;
                x[index + 2 * arraySize] = JMKCol;
            }
        }
//        for (int i = 0; i < x.length; i++) {
        for (int i = 0; i < arraySize; i++) {
            x[i] /= nPairs[i];
        }
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public CorrelationAgent makeAgent(IAtom a, Box agentBox) {

        return new CorrelationAgent(space);
    }

    public void releaseAgent(CorrelationAgent agent, IAtom atom, Box agentBox) {

    }

    public static class CorrelationAgent implements Integrator.Torquable, Integrator.Forcible {  //need public so to use with instanceof
        public final Vector torque;
        public final Vector force;

        public CorrelationAgent(Space space) {
            torque = new Vector1D();
            force = space.makeVector();
        }

        public Vector torque() {
            return torque;
        }

        public Vector force() {
            return force;
        }
    }
}
