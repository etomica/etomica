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
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.BondingInfo;
import etomica.potential.P2SoftSphere;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMasterFasterer;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;

/**
 * Simple soft-sphere Monte Carlo simulation in 3D.
 *
 * @author Tai Boon Tan
 */
public class SoftSphere3dFasterer extends Simulation {

    public IntegratorMCFasterer integrator;
    public MCMoveAtomFasterer mcMoveAtom;
    public SpeciesGeneral species;
    public Box box;
    public P2SoftSphere potential;
    public PotentialMasterFasterer potentialMaster;
    public DataSourceCountSteps meterCycles;

    public SoftSphere3dFasterer() {
        this(1.338, 12, .1);
    }

    public SoftSphere3dFasterer(double density, int exponent, double temperature) {
        super(Space3D.getInstance());

        species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(this));
        addSpecies(species);

        box = this.makeBox();
        potentialMaster = new PotentialMasterFasterer(getSpeciesManager(), box, BondingInfo.noBonding());
        integrator = new IntegratorMCFasterer(this, potentialMaster, box);
        integrator.setTemperature(temperature);


        mcMoveAtom = new MCMoveAtomFasterer(random, potentialMaster, box);

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
        potentialMaster.setPairPotential(type1, type1, truncated);
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

        final SoftSphere3dFasterer sim = new SoftSphere3dFasterer(density, exponent, temperature);
        int numAtoms = sim.box.getNMolecules(sim.species);

        MeterPotentialEnergyFromIntegratorFasterer meterEnergy = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
        AccumulatorAverageCollapsing accumulator = new AccumulatorAverageCollapsing();
        DataPumpListener pump = new DataPumpListener(meterEnergy, accumulator);

        accumulator.setPushInterval(1);
        sim.integrator.getEventManager().addListener(pump);

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, 10000000));


        double temp = sim.integrator.getTemperature();
        double Cv = ((DataDouble) ((DataGroup) accumulator.getData()).getData(accumulator.STANDARD_DEVIATION.index)).x;
        double energy = ((DataDouble) ((DataGroup) accumulator.getData()).getData(accumulator.AVERAGE.index)).x;
        Cv /= temp;
        Cv *= Cv / numAtoms;
        System.out.println("Cv/k: " + Cv);
        System.out.println("System Energy: " + energy / numAtoms);

    }


}
