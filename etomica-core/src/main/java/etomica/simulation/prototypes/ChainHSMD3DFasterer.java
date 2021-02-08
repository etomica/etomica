/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.simulation.prototypes;


import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomPair;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.config.ConformationLinear;
import etomica.data.AccumulatorHistogram;
import etomica.data.DataPumpListener;
import etomica.data.histogram.HistogramCollapsing;
import etomica.data.meter.MeterRadiusGyration;
import etomica.graphics.ColorSchemeRandomByMolecule;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHardFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.IMolecule;
import etomica.nbr.list.NeighborListManagerFastererHard;
import etomica.potential.P2HardGeneric;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialComputePair;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesBuilder;
import etomica.species.SpeciesGeneral;

import java.util.ArrayList;
import java.util.List;

/**
 * Molecular dynamics of chains of hard spheres.
 * <p>
 * Each chain has 4 spheres, joined by a hard tether potential.
 */
public class ChainHSMD3DFasterer extends Simulation {

    public Box box;
    public IntegratorHardFasterer integrator;
    public SpeciesGeneral species;
    public P2HardGeneric potential;
    public AccumulatorHistogram histogramRG;
    public DataPumpListener pumpRG;

    public ChainHSMD3DFasterer() {
        super(Space3D.getInstance());

        int chainLength = 4;
        double sigma = 1.0;
        double bondFactor = 0.15;
        P2HardGeneric potentialBond = new P2HardGeneric(new double[]{sigma - bondFactor, sigma + bondFactor, 2 * (sigma + bondFactor)},
                new double[]{Double.POSITIVE_INFINITY, 0, Double.POSITIVE_INFINITY});
        ConformationLinear conformation = new ConformationLinear(space);
        conformation.setBondLength(1.0);
        conformation.setAngle(1, 0.35);
        species = new SpeciesBuilder(space)
                .withConformation(conformation)
                .setDynamic(true)
                .addCount(AtomType.simpleFromSim(this), chainLength)
                .build();
        this.addSpecies(species);

        double neighborRangeFac = 1.6;
        PotentialMasterBonding.FullBondingInfo bondingInfo = new PotentialMasterBonding.FullBondingInfo(this);
        box = this.makeBox();
        NeighborListManagerFastererHard neighborManager = new NeighborListManagerFastererHard(getSpeciesManager(), box, 2, neighborRangeFac, bondingInfo);
        neighborManager.setDoDownNeighbors(true);
        PotentialComputePair potentialMaster = new PotentialComputePair(this, box, neighborManager);

        int numAtoms = 108;

        integrator = new IntegratorHardFasterer(potentialMaster, neighborManager, random, 0.01, 1.0, box, bondingInfo);
        integrator.setIsothermal(false);

        getController().addActivity(new ActivityIntegrate(integrator, true));

        List<int[]> bondedIndices = new ArrayList<>();
        for (int i = 0; i < chainLength - 1; i++) {
            bondedIndices.add(new int[]{i, i + 1});
        }
        bondingInfo.setBondingPotentialPair(species, potentialBond, bondedIndices);

        double l = 14.4573 * Math.pow((chainLength * numAtoms / 2020.0), 1.0 / 3.0);
        box.getBoundary().setBoxSize(Vector.of(new double[]{l, l, l}));
        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        box.setNMolecules(species, numAtoms);
        config.initializeCoordinates(box);

        this.potential = P2HardSphere.makePotential(sigma);
        AtomType leafType = species.getAtomType(0);
        potentialMaster.setPairPotential(leafType, leafType, potential);

        MeterRadiusGyration meterRG = new MeterRadiusGyration(space);
        meterRG.setBox(box);
        histogramRG = new AccumulatorHistogram(new HistogramCollapsing(), 10);
        pumpRG = new DataPumpListener(meterRG, histogramRG);
        integrator.getEventManager().addListener(pumpRG);
    }

    public static void main(String[] args) {

        final ChainHSMD3DFasterer sim = new ChainHSMD3DFasterer();
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
