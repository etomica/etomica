package etomica.spin.heisenberg;

import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPairSum extends DataSourceScalar {


    protected final Space space;
    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective allAtoms;
    protected double temperature;
    protected double J;
    protected double mu;
    protected double bt;
    protected Vector dr;
    protected Vector work;
    protected AtomLeafAgentManager leafAgentManager;
    private Box box;

    //TODO debug only
    protected PotentialCalculationSumquare Ans;

    public MeterPairSum(Box box, double temperature, double interactionS, double dipoleMagnitude, PotentialMaster potentialMaster) {
        super("anything", Null.DIMENSION);
        this.box = box;
        this.space = box.getSpace();
        this.temperature = temperature;
        this.potentialMaster = potentialMaster;
        J = interactionS;
        bt = 1 / temperature;
        mu = dipoleMagnitude;


        Ans = new PotentialCalculationSumquare(space, dipoleMagnitude);
        allAtoms = new IteratorDirective();

    }

    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("no box");
        Ans.zeroSum();
        potentialMaster.calculate(box, allAtoms, Ans);

        return Ans.getSum();
    }

}
