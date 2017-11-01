package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.LatticeCubicFcc;
import etomica.listener.IntegratorListenerAction;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

/**
 * A demo simulation to test aspects of the frontend.
 * <p>
 * Separate paragraph
 * <p>
 * Another paragraph. {@code some.code;}
 * <p>
 * <ul>
 * <li>List item</li>
 * </ul>
 */
public class DemoSim extends Simulation {

    private final PotentialMasterMonatomic potentialMaster;

    public DemoSim() {
        super(Space3D.getInstance());
        potentialMaster = new PotentialMasterMonatomic(this);
        IntegratorMC integrator = new IntegratorMC(this, potentialMaster);
        integrator.setTemperature(1);


        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(10000000);
        getController().addAction(activityIntegrate);

        SpeciesSpheresMono species = new SpeciesSpheresMono(this, space);
        //species2 = new SpeciesSpheresMono(this);
        addSpecies(species);
        Box box = new Box(space);
        addBox(box);
        box.setNMolecules(species, 108);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(1);
        inflater.actionPerformed();
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        P2SoftSphere potential = new P2SoftSphere(space, 1, 1, 12);
        P2SoftSphericalTruncated truncated = new P2SoftSphericalTruncated(space, potential, box.getBoundary().getBoxSize().getX(0) / 2);

        AtomType type1 = species.getLeafType();
        potentialMaster.addPotential(truncated, new AtomType[]{type1, type1});

        DataSourceCountSteps meterCycles = new DataSourceCountSteps(integrator);

        integrator.setBox(box);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));

        MCMoveVolume mcMoveVolume = new MCMoveVolume(potentialMaster, this.getRandom(), this.getSpace(), 1);
        integrator.getMoveManager().addMCMove(mcMoveVolume);

        MeterPotentialEnergy energy = new MeterPotentialEnergy(potentialMaster);
        energy.setBox(box);

        AccumulatorHistogram hist = new AccumulatorHistogram();
        DataPump histPump = new DataPump(energy, hist);
        IntegratorListenerAction histAct = new IntegratorListenerAction(histPump);
        integrator.getEventManager().addListener(histAct);

        AccumulatorAverageCollapsing avgEnergy = new AccumulatorAverageCollapsing();
        avgEnergy.setPushInterval(10);
        DataPump pump = new DataPump(energy, avgEnergy);
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(pump);
        pumpListener.setInterval(10);
        integrator.getEventManager().addListener(pumpListener);

        AccumulatorHistory h = new AccumulatorHistory(new HistoryCollapsingAverage());
        DataPumpListener l = new DataPumpListener(energy, h, 100);
        integrator.getEventManager().addListener(l);

    }

}
