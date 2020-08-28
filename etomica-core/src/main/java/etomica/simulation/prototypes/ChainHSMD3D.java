/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.iterator.ApiBuilder;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.meter.MeterRadiusGyration;
import etomica.graphics.*;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.nbr.CriterionAll;
import etomica.nbr.CriterionBondedSimple;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardBond;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

/**
 * Molecular dynamics of chains of hard spheres.
 * <p>
 * Each chain has 4 spheres, joined by a hard tether potential.
 */
public class ChainHSMD3D extends Simulation {

    public Box box;
    public IntegratorHard integrator;
    public SpeciesGeneral species;
    public P2HardSphere potential;
    public AccumulatorHistogram histogramRG;
    public DataPumpListener pumpRG;

    public ChainHSMD3D() {
        super(Space3D.getInstance());

        int chainLength = 4;
        P2HardBond potentialBond = new P2HardBond(space, 1.0, 0.15, true);
        ConformationLinear conformation = new ConformationLinear(space);
        conformation.setBondLength(1.0);
        conformation.setAngle(1, 0.35);
        species = new SpeciesBuilder(space)
                .withConformation(conformation)
                .setDynamic(true)
                .addCount(AtomType.simpleFromSim(this), chainLength)
                .build();
        this.addSpecies(species);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, space);
        int numAtoms = 108;
        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac);

        box = this.makeBox();
        integrator = new IntegratorHard(this, potentialMaster, box);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        getController().addActivity(new ActivityIntegrate(integrator, true));

        PotentialGroup group = potentialMaster.makePotentialGroup(1);
        group.addPotential(potentialBond, ApiBuilder.makeAdjacentPairIterator());
        potentialMaster.addPotential(group, new ISpecies[]{species});

        double l = 14.4573 * Math.pow((chainLength * numAtoms / 2020.0), 1.0 / 3.0);
        box.getBoundary().setBoxSize(Vector.of(new double[]{l, l, l}));
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        box.setNMolecules(species, numAtoms);
        config.initializeCoordinates(box);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        this.potential = new P2HardSphere(space, 1.0, true);
        AtomType leafType = species.getAtomType(0);
        potentialMaster.addPotential(this.potential, new AtomType[]{leafType, leafType});
        CriterionBondedSimple nonBondedCriterion = new CriterionBondedSimple(new CriterionAll());
        nonBondedCriterion.setBonded(false);
        ((CriterionInterMolecular) potentialMaster.getCriterion(leafType, leafType)).setIntraMolecularCriterion(nonBondedCriterion);

        MeterRadiusGyration meterRG = new MeterRadiusGyration(space);
        meterRG.setBox(box);
        histogramRG = new AccumulatorHistogram(new HistogramCollapsing(), 10);
        pumpRG = new DataPumpListener(meterRG, histogramRG);
        integrator.getEventManager().addListener(pumpRG);
    }

    public static void main(String[] args) {

        final etomica.simulation.prototypes.ChainHSMD3D sim = new etomica.simulation.prototypes.ChainHSMD3D();
        final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
        DisplayPlot plotRG = new DisplayPlot();
        sim.histogramRG.addDataSink(plotRG.getDataSet().makeDataSink());
        plotRG.setLabel("RG");
        plotRG.setDoLegend(false);
        simGraphic.add(plotRG);
        simGraphic.getController().getDataStreamPumps().add(sim.pumpRG);
        for (IMolecule molecule : sim.box.getMoleculeList(sim.species)) {
            for (int i = 1; i < molecule.getChildList().size(); i++) {
                IAtom atom1 = molecule.getChildList().get(i - 1);
                IAtom atom2 = molecule.getChildList().get(i);
                ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).makeBond(new AtomPair(atom1, atom2), null);
            }
        }
        ColorSchemeRandomByMolecule colorScheme = new ColorSchemeRandomByMolecule(sim, sim.box, sim.getRandom());
        simGraphic.getDisplayBox(sim.box).setColorScheme(colorScheme);

        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        simGraphic.makeAndDisplayFrame();
    }
}
