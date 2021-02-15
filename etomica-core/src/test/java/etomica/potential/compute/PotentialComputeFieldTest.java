package etomica.potential.compute;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceScalar;
import etomica.data.meter.MeterKineticEnergy;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.IntegratorVelocityVerletFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.potential.IPotentialField;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.units.dimensions.Null;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class PotentialComputeFieldTest {

    Simulation sim;
    IntegratorMCFasterer integrator;
    IntegratorVelocityVerletFasterer integratorMD;
    AccumulatorAverageFixed acc, accKE;
    DataPumpListener pump;

    @BeforeEach
    void setUp() {
        sim = new Simulation(Space3D.getInstance());
        ISpecies species = SpeciesGeneral.monatomic(sim.getSpace(), AtomType.simple("A"), true);
        sim.addSpecies(species);
        Box box = sim.makeBox();
        box.getBoundary().setBoxSize(Vector.of(10, 10, 10));
        box.setNMolecules(species, 100);
        PotentialComputeField pcOne = new PotentialComputeField(sim.getSpeciesManager(), box);
        IPotentialField p1 = new IPotentialField() {
            @Override
            public double u(IAtom atom) {
                double x = atom.getPosition().getX(0);
                return 0.5 * x * x;
            }

            @Override
            public double udu(IAtom atom, Vector f) {
                double x = atom.getPosition().getX(0);
                if (Double.isNaN(x)) throw new RuntimeException("oops");
                double fx = f.getX(0);
                f.setX(0, fx - x);
                if (f.isNaN()) throw new RuntimeException("oops");
                return 0.5 * x * x;
            }
        };
        pcOne.setFieldPotential(species.getLeafType(), p1);
        integratorMD = new IntegratorVelocityVerletFasterer(pcOne, sim.getRandom(), 0.0005, 1.0, box);
        integratorMD.setIsothermal(true);
        integratorMD.setThermostatInterval(10);
        integrator = new IntegratorMCFasterer(pcOne, sim.getRandom(), 1, box);
        integrator.getMoveManager().addMCMove(new MCMoveAtomFasterer(sim.getRandom(), pcOne, box));
        sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, 10000));

        DataSourceScalar msx = new DataSourceScalar("msx", Null.DIMENSION) {
            @Override
            public double getDataAsScalar() {
                double sum = 0;
                IAtomList atoms = box.getLeafList();
                for (IAtom a : atoms) {
                    double x = a.getPosition().getX(0);
                    if (Double.isNaN(x)) throw new RuntimeException("oops");
                    sum += x * x;
                }
                return sum / atoms.size();
            }
        };
        acc = new AccumulatorAverageFixed(1);
        pump = new DataPumpListener(msx, acc, 100);
        integrator.getEventManager().addListener(pump);
        integratorMD.getEventManager().addListener(pump);

        MeterKineticEnergy meterKE = new MeterKineticEnergy(box);
        accKE = new AccumulatorAverageFixed(1);
        DataPumpListener pumpKE = new DataPumpListener(meterKE, accKE, 10);
        integratorMD.getEventManager().addListener(pumpKE);
    }

    @Test
    public void testStdevMC() {
        acc.reset();
        pump.setInterval(100);
        sim.getController().runActivityBlocking(new ActivityIntegrate(integrator, 100000));
        double stdev = acc.getData(acc.AVERAGE).getValue(0);
        Assertions.assertEquals(stdev, 1, 0.05);
    }

    @Test
    public void testStdevMD() {
        pump.setInterval(10);
        acc.reset();
        sim.getController().runActivityBlocking(new ActivityIntegrate(integratorMD, 100000));
        double stdev = acc.getData(acc.AVERAGE).getValue(0);
        // do we really do 100x more work and still need 8x stdev?  yes!
        Assertions.assertEquals(stdev, 1, 0.4);
    }
}