/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.DataPump;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

/**
 * Simple soft-sphere Monte Carlo simulation in 3D.
 *
 * @author Tai Boon Tan
 */
public class SoftSphere3d extends Simulation {

    public IntegratorMC integrator;
    public MCMoveAtom mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public P2SoftSphere potential;
    public PotentialMaster potentialMaster;
    public DataSourceCountSteps meterCycles;

    public SoftSphere3d() {
        this(1.338, 12, .1);
    }

    public SoftSphere3d(double density, int exponent, double temperature) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        potentialMaster = new PotentialMasterMonatomic(getSpeciesManager());
        box = this.makeBox();
        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setTemperature(temperature);


        mcMoveAtom = new MCMoveAtom(random, potentialMaster, space);
        getController().addActivity(new ActivityIntegrate(integrator), 10000000);

        box.setNMolecules(species, 108);
        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(density);
        inflater.actionPerformed();
        // box.setNMolecules(species2, 20);
        new ConfigurationLattice(new LatticeCubicFcc(space), space).initializeCoordinates(box);
        potential = new P2SoftSphere(space, 1, 1, exponent);
        P2SoftSphericalTruncated truncated = new P2SoftSphericalTruncated(space, potential, box.getBoundary().getBoxSize().getX(0) / 2);
        // System.out.println("Truncated radius is: " +truncated.getTruncationRadius());

        AtomType type1 = species.getLeafType();
        //AtomType type2 = species2.getLeafType();
        potentialMaster.addPotential(truncated, new AtomType[]{type1, type1});
        // potentialMaster.addPotential(potential, new AtomType[] {type1, type2});
        //potentialMaster.addPotential(potential, new AtomType[] {type2, type2});

        meterCycles = new DataSourceCountSteps(integrator);
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getEventManager().addListener(new IntegratorListenerAction(new BoxImposePbc(box, space)));

//	    LatticeRenderer.ColorSchemeCell colorSchemeCell = new LatticeRenderer.ColorSchemeCell();
//	    display.setColorScheme(colorSchemeCell);

//		elementCoordinator.go();
//	    etomica.lattice.BravaisLattice lattice = ((IteratorFactoryCell)this.getIteratorFactory()).getLattice(box);
//        colorSchemeCell.setLattice(lattice);
    }

    public static void main(String[] args) {

        double density = 1.338;
        int exponent = 12;
        double temperature = 0.1;


        if (args.length > 0) {
            density = Double.parseDouble(args[0]);
        }
        if (args.length > 1) {
            exponent = Integer.parseInt(args[1]);
        }
        if (args.length > 1) {
            temperature = Double.parseDouble(args[2]);
        }

        final SoftSphere3d sim = new SoftSphere3d(density, exponent, temperature);
        int numAtoms = sim.box.getNMolecules(sim.species);

        MeterPotentialEnergyFromIntegrator meterEnergy = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        DataPump pump = new DataPump(meterEnergy, null);
        AccumulatorAverageCollapsing accumulator = new AccumulatorAverageCollapsing();

        accumulator.setPushInterval(1);
        pump.setDataSink(accumulator);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(pump));

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, Long.MAX_VALUE));


        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) accumulator.getData()).getData(accumulator.STANDARD_DEVIATION.index)).x;
        double energy = ((DataDouble) ((DataGroup) accumulator.getData()).getData(accumulator.AVERAGE.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k: " + Cv);
        System.out.println("System Energy: " + energy / numAtoms);

    }


}
