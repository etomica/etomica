/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.paracetamol;

import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.normalmode.MCMoveHarmonicStep;
import etomica.normalmode.NormalModesFromFile;
import etomica.normalmode.WaveVectorFactory;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 * 
 * Orthorhombic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolOrthorhombicMoveHarmonicStep extends Simulation {

	private static final long serialVersionUID = 1L;
	//private final static String APP_NAME = "MC Move Harmonic Step Paracetamol Orthorhombic";
    public Box box;
    public IntegratorMC integrator;
    public SpeciesParacetamol species;
    public P2ElectrostaticDreiding potentialCC , potentialCHy , potentialHyHy;
    public P2ElectrostaticDreiding potentialCN , potentialNO  , potentialNN  ;
    public P2ElectrostaticDreiding potentialHyN, potentialHyO , potentialOO  ;
    public P2ElectrostaticDreiding potentialCO , potentialHpHp, potentialCHp ;
    public P2ElectrostaticDreiding potentialHpN, potentialOHp , potentialHyHp;
    public Controller controller;
    public PotentialMaster potentialMaster;
    public BravaisLattice lattice;
    public BoundaryRectangularPeriodic bdry;
    public ActivityIntegrate actionIntegrate;

    public MCParacetamolOrthorhombicMoveHarmonicStep() {
        this(192);
    }

    private MCParacetamolOrthorhombicMoveHarmonicStep(int numMolecules) {
        super(Space3D.getInstance());

        species = new SpeciesParacetamol(space, false);
        addSpecies(species);

        potentialMaster = new PotentialMaster();

        /*
         * Orthorhombic Crystal
         */

        PrimitiveOrthorhombic primitive = new PrimitiveOrthorhombic(space, 17.248, 12.086, 7.382);
        // 17.248, 12.086, 7.382
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();
        lattice = new BravaisLatticeCrystal(primitive, basis);

        bdry = new BoundaryRectangularPeriodic(space, 1); //unit cell
        bdry.setBoxSize(space.makeVector(new double[]{2 * 17.248, 3 * 12.086, 4 * 7.382}));
        box = this.makeBox(bdry);
        integrator = new IntegratorMC(this, potentialMaster, box);
        integrator.setIsothermal(false);
        integrator.setTemperature(Kelvin.UNIT.toSim(1));


        /*
         * Harmonic Simulation
         */
        MCMoveHarmonicStep moveHarmonicStep = new MCMoveHarmonicStep(potentialMaster, random);
        integrator.getMoveManager().addMCMove(moveHarmonicStep);

        NormalModesFromFile normalModes = new NormalModesFromFile("Normal_Modes_Paracetamol_FormII_100.0K", 3);
        normalModes.setTemperature(Kelvin.UNIT.toSim(100));

        ConformationParacetamolOrthorhombic conformation = new ConformationParacetamolOrthorhombic(space);
        species.setConformation(conformation);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{25, 25, 25}));
        box.setNMolecules(species, numMolecules);

        WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
        waveVectorFactory.makeWaveVectors(box);

        moveHarmonicStep.setBox(box);
        int[] num = new int[45];
        for (int i = 0; i < 45; i++) {
            num[i] = i + 3;
        }
        moveHarmonicStep.setModes(num);
        moveHarmonicStep.setEigenVectors(normalModes.getEigenvectors()[0]);
        moveHarmonicStep.setStepSize(0.01);
        ((MCMoveStepTracker) moveHarmonicStep.getTracker()).setAdjustInterval(2);
        ((MCMoveStepTracker) moveHarmonicStep.getTracker()).setNoisyAdjustment(true);


        actionIntegrate = new ActivityIntegrate(integrator, 0, false);
        actionIntegrate.setMaxSteps(1000000);
        getController().addAction(actionIntegrate);

        /*
         * Intermolecular Potential
         */
        double truncationRadiusCC = 3.0 * 3.395524116;
        double truncationRadiusCHy = 3.0 * 2.670105986;
        double truncationRadiusHyHy = 3.0 * 2.099665865;
        double truncationRadiusCN = 3.0 * 3.237739512;
        double truncationRadiusNO = 3.0 * 3.035146951;
        double truncationRadiusNN = 3.0 * 3.087286897;
        double truncationRadiusHyN = 3.0 * 2.546030404;
        double truncationRadiusHyO = 3.0 * 2.503031484;
        double truncationRadiusOO = 3.0 * 2.983887553;
        double truncationRadiusCO = 3.0 * 3.183058614;
        double truncationRadiusHpHp = 3.0 * 1.543178334;
        double truncationRadiusCHp = 3.0 * 2.289082817;
        double truncationRadiusHpN = 3.0 * 2.182712943;
        double truncationRadiusOHp = 3.0 * 2.145849932;
        double truncationRadiusHyHp = 3.0 * 1.800044389;

        potentialCC = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3832.14700),
                0.277778, ElectronVolt.UNIT.toSim(25.286949));
        potentialCHy = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(689.53672),
                0.272480, ElectronVolt.UNIT.toSim(5.978972));
        potentialHyHy = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(124.07167),
                0.267380, ElectronVolt.UNIT.toSim(1.413698));
        potentialCN = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3179.51460),
                0.271003, ElectronVolt.UNIT.toSim(19.006710));
        potentialNO = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2508.04480),
                0.258398, ElectronVolt.UNIT.toSim(12.898341));
        potentialNN = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2638.02850),
                0.264550, ElectronVolt.UNIT.toSim(14.286224));
        potentialHyN = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(572.10541),
                0.265957, ElectronVolt.UNIT.toSim(4.494041));
        potentialHyO = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(543.91604),
                0.259740, ElectronVolt.UNIT.toSim(4.057452));
        potentialOO = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2384.46580),
                0.252525, ElectronVolt.UNIT.toSim(11.645288));
        potentialCO = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3022.85020),
                0.264550, ElectronVolt.UNIT.toSim(17.160239));
        potentialHpHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(52.12899),
                0.214592, ElectronVolt.UNIT.toSim(0.222819));
        potentialCHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(446.95185),
                0.242131, ElectronVolt.UNIT.toSim(2.373693));
        potentialHpN = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(370.83387),
                0.236967, ElectronVolt.UNIT.toSim(1.784166));
        potentialOHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(352.56176),
                0.232019, ElectronVolt.UNIT.toSim(1.610837));
        potentialHyHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(80.42221),
                0.238095, ElectronVolt.UNIT.toSim(0.561248));

        // CA-CA
        if (truncationRadiusCC > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialCC = new P2SoftSphericalTruncated(space, potentialCC, truncationRadiusCC);
        potentialMaster.addPotential(interpotentialCC, new AtomType[]{species.getCType(), species.getCType()});

        // CA-HY
        if (truncationRadiusCHy > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialCHy = new P2SoftSphericalTruncated(space, potentialCHy, truncationRadiusCHy);
        potentialMaster.addPotential(interpotentialCHy, new AtomType[]{species.getCType(), species.getHyType()});

        // HY-HY
        if (truncationRadiusHyHy > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialHyHy = new P2SoftSphericalTruncated(space, potentialHyHy, truncationRadiusHyHy);
        potentialMaster.addPotential(interpotentialHyHy, new AtomType[]{species.getHyType(), species.getHyType()});

        // CA-NI
        if (truncationRadiusCN > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialCN = new P2SoftSphericalTruncated(space, potentialCN, truncationRadiusCN);
        potentialMaster.addPotential(interpotentialCN, new AtomType[]{species.getCType(), species.getNType()});

        // NI-OX
        if (truncationRadiusNO > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialNO = new P2SoftSphericalTruncated(space, potentialNO, truncationRadiusNO);
        potentialMaster.addPotential(interpotentialNO, new AtomType[]{species.getNType(), species.getOType()});

        //NI-NI
        if (truncationRadiusNN > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialNN = new P2SoftSphericalTruncated(space, potentialNN, truncationRadiusNN);
        potentialMaster.addPotential(interpotentialNN, new AtomType[]{species.getNType(), species.getNType()});

        // HY-NI
        if (truncationRadiusHyN > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialHyN = new P2SoftSphericalTruncated(space, potentialHyN, truncationRadiusHyN);
        potentialMaster.addPotential(interpotentialHyN, new AtomType[]{species.getHyType(), species.getNType()});

        // HY-OX
        if (truncationRadiusHyO > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large. " +
                    " Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialHyO = new P2SoftSphericalTruncated(space, potentialHyO, truncationRadiusHyO);
        potentialMaster.addPotential(interpotentialHyO, new AtomType[]{species.getHyType(), species.getOType()});

        // OX-OX
        if (truncationRadiusOO > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large. " +
                    " Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialOO = new P2SoftSphericalTruncated(space, potentialOO, truncationRadiusOO);
        potentialMaster.addPotential(interpotentialOO, new AtomType[]{species.getOType(), species.getOType()});

        // CA-OX
        if (truncationRadiusCO > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large. " +
                    " Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialCO = new P2SoftSphericalTruncated(space, potentialCO, truncationRadiusCO);
        potentialMaster.addPotential(interpotentialCO, new AtomType[]{species.getCType(), species.getOType()});

        // HP-HP
        if (truncationRadiusHpHp > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large. " +
                    " Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialHpHp = new P2SoftSphericalTruncated(space, potentialHpHp, truncationRadiusHpHp);
        potentialMaster.addPotential(interpotentialHpHp, new AtomType[]{species.getHpType(), species.getHpType()});

        // CA-HP
        if (truncationRadiusCHp > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large. " +
                    " Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialCHp = new P2SoftSphericalTruncated(space, potentialCHp, truncationRadiusCHp);
        potentialMaster.addPotential(interpotentialCHp, new AtomType[]{species.getCType(), species.getHpType()});

        // HP-NI
        if (truncationRadiusHpN > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialHpN = new P2SoftSphericalTruncated(space, potentialHpN, truncationRadiusHpN);
        potentialMaster.addPotential(interpotentialHpN, new AtomType[]{species.getHpType(), species.getNType()});

        // OX-HP
        if (truncationRadiusOHp > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large. " +
                    " Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialOHp = new P2SoftSphericalTruncated(space, potentialOHp, truncationRadiusOHp);
        potentialMaster.addPotential(interpotentialOHp, new AtomType[]{species.getOType(), species.getHpType()});

        // HY-HP
        if (truncationRadiusHyHp > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
                    "Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        P2SoftSphericalTruncated interpotentialHyHp = new P2SoftSphericalTruncated(space, potentialHyHp, truncationRadiusHyHp);
        potentialMaster.addPotential(interpotentialHyHp, new AtomType[]{species.getHyType(), species.getHpType()});

        potentialMaster.lrcMaster().setEnabled(false);
        /*
         *
         */


        CoordinateDefinitionParacetamol coordDef = new CoordinateDefinitionParacetamol(this, box, primitive, basis, space);
        coordDef.setBasisOrthorhombic();
        coordDef.initializeCoordinates(new int[]{2, 3, 4});

        moveHarmonicStep.setCoordinateDefinition(coordDef);


    } //end of constructor

    public static void main(String[] args) {
    	int numMolecules = 192;
        MCParacetamolOrthorhombicMoveHarmonicStep sim = new MCParacetamolOrthorhombicMoveHarmonicStep(numMolecules);
//        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME, 1);
//        Pixel pixel = new Pixel(10);
//        simGraphic.getDisplayPhase(sim.phase).setPixelUnit(pixel);

//        ConfigurationFile configFile = new ConfigurationFile("Coord_Paracetamol_FormII_100.0_K");
//        configFile.initializeCoordinates(sim.phase);
        /*****************************************************************************/

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
//        DisplayBox PEbox = new DisplayBox();
//        DataPump PEpump = new DataPump(meterPE, PEbox);

        WriteConfiguration writeConfig = new WriteConfiguration(sim.space);
        writeConfig.setConfName("Coord_Paracetamol_FormII_Minimum_Energy");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);

        IntegratorListenerAction writeConfigListener = new IntegratorListenerAction(writeConfig);
        writeConfigListener.setInterval(100);
        sim.integrator.getEventManager().addListener(writeConfigListener);


        /**********************************************************************/
//        simGraphic.add(PEbox);

        sim.getController().actionPerformed();

//        simGraphic.makeAndDisplayFrame(APP_NAME);
//        simGraphic.getDisplayPhase(sim.phase).setPixelUnit(new Pixel(10));
//        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.
//        		displayList().getFirst()).getColorScheme());
//        AtomTypeGroup atomType = (AtomTypeGroup)sim.species.getMoleculeType();
//        colorScheme.setColor(atomType.getChildTypes()[0], java.awt.Color.red);
//        colorScheme.setColor(atomType.getChildTypes()[1], java.awt.Color.gray);
//        colorScheme.setColor(atomType.getChildTypes()[2], java.awt.Color.blue);
//        colorScheme.setColor(atomType.getChildTypes()[3], java.awt.Color.white);
//        colorScheme.setColor(atomType.getChildTypes()[4], java.awt.Color.white);
//
//        simGraphic.getDisplayPhase(sim.phase).repaint();


    }//end of main
}//end of class
