package etomica.spin.heisenberg;

import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterMappedAveragingFreeEnergy extends DataSourceScalar {
    protected final PotentialCompute potentialMaster;
    protected final NeighborIterator nbrIterator;
    protected double J;
    protected double mu;
    protected double bt;
    private final Box box;

    public MeterMappedAveragingFreeEnergy(Box box, double temperature, double interactionS, double dipoleMagnitude, PotentialCompute potentialMaster, NeighborManager nbrManager) {
        super("Stuff", Null.DIMENSION);
        this.box = box;
        this.potentialMaster = potentialMaster;
        this.nbrIterator = nbrManager.makeNeighborIterator();
        J = interactionS;
        bt = 1 / temperature;
        mu = dipoleMagnitude;
    }

    public double getDataAsScalar() {
        double[] x = new double[]{0};
        potentialMaster.computeAll(true);
        Vector[] torques = potentialMaster.getTorques();
        for (IAtom a : box.getLeafList()) {
            Vector ei = ((IAtomOriented) a).getOrientation().getDirection();
            nbrIterator.iterUpNeighbors(a.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    IAtomOriented atom2 = (IAtomOriented) jAtom;
                    Vector ej = atom2.getOrientation().getDirection();


                    double cost1 = ei.getX(0);
                    double sint1 = ei.getX(1);
                    double sint2 = ej.getX(1);
                    double cost2 = ej.getX(0);
                    double sint1mt2 = sint1 * cost2 - cost1 * sint2;

                    double f1 = bt * torques[a.getLeafIndex()].getX(0);
                    double f2 = bt * torques[jAtom.getLeafIndex()].getX(0);
                    x[0] += 0.5 * J * (f1 - f2) * sint1mt2;
                }
            });
        }
        return x[0];
    }
}
