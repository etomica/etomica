/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


package etomica.paracetamol;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.iterator.ApiIndexList;
import etomica.atom.iterator.Atomset3IteratorIndexList;
import etomica.atom.iterator.Atomset4IteratorIndexList;
import etomica.box.Box;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.integrator.IntegratorListenerAction;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.CriterionNone;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.units.Kelvin;
import etomica.units.Pixel;

import java.awt.*;
import java.util.ArrayList;

/**
 * 
 * Three-dimensional soft-sphere molecular dynamics simulation for paracetamol molecule, using
 * neighbor listing.  
 * 
 * Orthorhombic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MDParacetamolOrthorhombic extends Simulation {

	private static final String APP_NAME = "MD Paracetamol Orthorhombic";
	private static final int PIXEL_SIZE = 15;

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    public final PotentialMasterList potentialMaster;
    
    /**
     * The Box holding the atoms. 
     */
    public final Box box;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorVelocityVerlet integrator;
    /**
     * The single soft-sphere species.
     */
    public final SpeciesParacetamol species;
    public BravaisLattice lattice;
    public BoundaryRectangularPeriodic bdry;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ConfigurationOrthorhombicLattice configOrthoLattice;
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    /*
     * we use a second, private constructor to permit the space to
     * appear twice in the call to the superclass constructor; alternatively
     * we could have passed Space3D.getInstance() twice
     *
     */
    private MDParacetamolOrthorhombic() {

        /*
         * invoke the superclass constructor
         *	"true" is indicating to the superclass that this is a dynamic simulation
         * the PotentialMaster is selected such as to implement neighbor listing
         */
        super(Space3D.getInstance());

        potentialMaster = new PotentialMasterList(this, 1.6, space);

        /*
         * Orthorhombic Crystal
         */

        PrimitiveOrthorhombic primitive = new PrimitiveOrthorhombic(space, 17.248, 12.086, 7.382);
        // 17.248, 12.086, 7.382
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();
        lattice = new BravaisLatticeCrystal(primitive, basis);
        configOrthoLattice = new ConfigurationOrthorhombicLattice(lattice, space);

        double neighborRangeFac = 1.6;

        box = new Box(space);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTimeStep(0.001); //1 = pico sec
        integrator.setTemperature(Kelvin.UNIT.toSim(123));

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        //activityIntegrate.setMaxSteps(100000);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        addBox(box);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{25, 25, 25}));
        species = new SpeciesParacetamol(space, true);
        addSpecies(species);
        box.setNMolecules(species, 128);

        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        PotentialGroup intramolecularpotential = potentialMaster.makePotentialGroup(1);
        potentialMaster.addPotential(intramolecularpotential, new ISpecies[]{species});

        /*
         *  Bond Stretch Potential
         *
         *	Equalibrium Radius [unit Amstrom]; Pre-factor [unit Kelvin]
         */

        P2Dreiding potentialC4O1 = new P2Dreiding(space, 1.352, .17612581e6);
        intramolecularpotential.addPotential(potentialC4O1,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexO[0]}}));

        P2Dreiding potentialring1 = new P2Dreiding(space, 1.395, .264188715e6);
        intramolecularpotential.addPotential(potentialring1,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexC[5]},
                        {SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexC[1]},
                        {SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[3]},
                        {SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[4]},
                }));

        P2Dreiding potentialring2 = new P2Dreiding(space, 1.385, .264188715e6);
        intramolecularpotential.addPotential(potentialring2,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[2]},
                        {SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[5]}
                }));

        P2Dreiding potentialC7C8 = new P2Dreiding(space, 1.503, .17612581e6);
        intramolecularpotential.addPotential(potentialC7C8,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexC[7]}}));

        P2Dreiding potentialC1N1 = new P2Dreiding(space, 1.394, .17612581e6);
        intramolecularpotential.addPotential(potentialC1N1,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0]}}));

        P2Dreiding potentialC7N1 = new P2Dreiding(space, 1.366, .17612581e6);
        intramolecularpotential.addPotential(potentialC7N1,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexN[0]}}));

        P2Dreiding potentialC7O2 = new P2Dreiding(space, 1.226, .35225162e6);
        intramolecularpotential.addPotential(potentialC7O2,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexO[1]}}));

        /*
         * Bond Angle Potential
         *
         * Equilibrium Angle [unit radians]; Pre-factor [unit Kelvin]
         */

        P3BondAngleDreiding potential120 = new P3BondAngleDreiding(space, 2.094395102, .3354777333e5);
        intramolecularpotential.addPotential(potential120,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[4]},
                        {SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[5]},
                        {SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexC[0]},
                        {SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexC[1]},
                        {SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[2]},
                        {SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[3]}
                }));

        P3BondAngleDreiding potentialO1C4C5 = new P3BondAngleDreiding(space, 2.14675498, .3354777333e5);
        intramolecularpotential.addPotential(potentialO1C4C5,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexO[0], SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[4]}}));

        P3BondAngleDreiding potentialO1C4C3 = new P3BondAngleDreiding(space, 2.042035225, .3354777333e5);
        intramolecularpotential.addPotential(potentialO1C4C3,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexO[0], SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[2]}}));

        P3BondAngleDreiding potentialC6C1N1 = new P3BondAngleDreiding(space, 2.059488517, .3354777333e5);
        intramolecularpotential.addPotential(potentialC6C1N1,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0]}}));

        P3BondAngleDreiding potentialC2C1N1 = new P3BondAngleDreiding(space, 2.129301687, .3354777333e5);
        intramolecularpotential.addPotential(potentialC2C1N1,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0]}}));

        P3BondAngleDreiding potentialN1C7C8 = new P3BondAngleDreiding(space, 1.989675347, .3354777333e5);
        intramolecularpotential.addPotential(potentialN1C7C8,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexN[0], SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexC[7]}}));

        P3BondAngleDreiding potentialN1C7O2 = new P3BondAngleDreiding(space, 2.181661565, .3354777333e5);
        intramolecularpotential.addPotential(potentialN1C7O2,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexN[0], SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexO[1]}}));

        P3BondAngleDreiding potentialO2C7C8 = new P3BondAngleDreiding(space, 2.111848395, .3354777333e5);
        intramolecularpotential.addPotential(potentialO2C7C8,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexO[1], SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexC[7]}}));

        P3BondAngleDreiding potentialC1N1C7 = new P3BondAngleDreiding(space, 2.251474735, .2742552176e5);
        intramolecularpotential.addPotential(potentialC1N1C7,
                new Atomset3IteratorIndexList(new int[][]{{SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0], SpeciesParacetamol.indexC[6]}}));

        /*
         * Torsional Angle Potential
         *
         * Equilibrium Phi [unit radians]; Pre-factor [unit Kelvin]; Periodicity [integer]
         */

        P4TorsionDreiding potentialCRCR = new P4TorsionDreiding(space, 3.141592654, 6290.2075, 2);
        intramolecularpotential.addPotential(potentialCRCR,
                new Atomset4IteratorIndexList(new int[][]{{SpeciesParacetamol.indexO[0], SpeciesParacetamol.indexC[3],
                        SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[5]},
                        {SpeciesParacetamol.indexO[0], SpeciesParacetamol.indexC[3],
                                SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[1]},
                        {SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[4],
                                SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexC[0]},
                        {SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[2],
                                SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[0]},
                        {SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[5],
                                SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0]},
                        {SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[1],
                                SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0]}

                }));

        P4TorsionDreiding potentialCRN3 = new P4TorsionDreiding(space, 0.0, 251.6083, 6);
        intramolecularpotential.addPotential(potentialCRN3,
                new Atomset4IteratorIndexList(new int[][]{{SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexC[0],
                        SpeciesParacetamol.indexN[0], SpeciesParacetamol.indexC[6]},
                        {SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[0],
                                SpeciesParacetamol.indexN[0], SpeciesParacetamol.indexC[6]}
                }));


        P4TorsionDreiding potentialN3C3 = new P4TorsionDreiding(space, 3.141592654, 503.2166, 3);
        intramolecularpotential.addPotential(potentialN3C3,
                new Atomset4IteratorIndexList(new int[][]{{SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0],
                        SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexO[1]},
                        {SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexN[0],
                                SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexC[7]}
                }));


        /*
         * Intra-nonbonded Potential
         *
         * Equilibrium Radius [unit Amstrom]; Pre-factor [unit Kelvin]
         */

        P2Exp6 potentialCC = new P2Exp6(space, 3832.147000 * 11604.45728, 0.277778, 25.286949 * 11604.45728);
        intramolecularpotential.addPotential(potentialCC,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[6]},
                        {SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexC[7]},
                        {SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[6]},
                        {SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexC[7]},
                        {SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[6]},
                        {SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexC[7]},
                        {SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexC[7]},
                        {SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexC[7]}
                }));

        P2Exp6 potentialCO = new P2Exp6(space, 3022.850200 * 11604.45728, 0.264550, 17.160239 * 11604.45728);
        intramolecularpotential.addPotential(potentialCO,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[0], SpeciesParacetamol.indexO[0]},
                        {SpeciesParacetamol.indexC[6], SpeciesParacetamol.indexO[0]},
                        {SpeciesParacetamol.indexC[7], SpeciesParacetamol.indexO[0]},
                        {SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexO[1]},
                        {SpeciesParacetamol.indexC[4], SpeciesParacetamol.indexO[1]},
                        {SpeciesParacetamol.indexC[2], SpeciesParacetamol.indexO[1]},
                        {SpeciesParacetamol.indexC[5], SpeciesParacetamol.indexO[1]},
                        {SpeciesParacetamol.indexC[1], SpeciesParacetamol.indexO[1]}
                }));

        P2Exp6 potentialON = new P2Exp6(space, 2508.044800 * 11604.45728, 0.258398, 12.898341 * 11604.45728);
        intramolecularpotential.addPotential(potentialON,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexO[0], SpeciesParacetamol.indexN[0]}}));

        P2Exp6 potentialCN = new P2Exp6(space, 3179.514600 * 11604.45728, 0.271003, 19.006710 * 11604.45728);
        intramolecularpotential.addPotential(potentialCN,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexC[3], SpeciesParacetamol.indexN[0]}}));

        P2Exp6 potentialO1O2 = new P2Exp6(space, 2384.465800 * 11604.45728, 0.252525, 11.645288 * 11604.45728);
        intramolecularpotential.addPotential(potentialO1O2,
                new ApiIndexList(new int[][]{{SpeciesParacetamol.indexO[0], SpeciesParacetamol.indexO[1]}}));



        /*
         * Intermolecular Potential
         */

        double truncationRadiusCC = 3.0 * 3.472990473;
        double truncationRadiusCO = 3.0 * 3.253072125;
        double truncationRadiusCN = 3.0 * 3.369296217;
        double truncationRadiusON = 3.0 * 3.149377868;
        double truncationRadiusOO = 3.0 * 3.033153776;
        double truncationRadiusNN = 3.0 * 3.262560196;

        if (truncationRadiusCC > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2 * truncationRadiusCC);
        P2SoftSphericalTruncated interpotentialCC = new P2SoftSphericalTruncated(space, new P2ElectrostaticDreiding(space, 3832.14700 * 11604.45728, 0.277778, 25.286949 * 11604.45728), truncationRadiusCC);
        potentialMaster.addPotential(interpotentialCC, new AtomType[]{species.getCType(), species.getCType()});

        if (truncationRadiusCO > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2 * truncationRadiusCO);
        P2SoftSphericalTruncated interpotentialCO = new P2SoftSphericalTruncated(space, new P2ElectrostaticDreiding(space, 3022.850200 * 11604.45728, 0.264550, 17.160239 * 11604.45728), truncationRadiusCO);
        potentialMaster.addPotential(interpotentialCO, new AtomType[]{species.getCType(), species.getOType()});

        if (truncationRadiusCN > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2 * truncationRadiusCN);
        P2SoftSphericalTruncated interpotentialCN = new P2SoftSphericalTruncated(space, new P2ElectrostaticDreiding(space, 3179.514600 * 11604.45728, 0.271003, 19.006710 * 11604.45728), truncationRadiusCN);
        potentialMaster.addPotential(interpotentialCN, new AtomType[]{species.getCType(), species.getNType()});

        if (truncationRadiusON > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2 * truncationRadiusON);
        P2SoftSphericalTruncated interpotentialON = new P2SoftSphericalTruncated(space, new P2ElectrostaticDreiding(space, 2508.044800 * 11604.45728, 0.258398, 12.898341 * 11604.45728), truncationRadiusON);
        potentialMaster.addPotential(interpotentialON, new AtomType[]{species.getOType(), species.getNType()});

        if (truncationRadiusOO > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2 * truncationRadiusOO);
        P2SoftSphericalTruncated interpotentialOO = new P2SoftSphericalTruncated(space, new P2ElectrostaticDreiding(space, 2384.465800 * 11604.45728, 0.252525, 11.645288 * 11604.45728), truncationRadiusOO);
        potentialMaster.addPotential(interpotentialOO, new AtomType[]{species.getOType(), species.getOType()});

        if (truncationRadiusNN > 0.5 * box.getBoundary().getBoxSize().getX(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is" + 0.5 * box.getBoundary().getBoxSize().getX(0));
        }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2 * truncationRadiusNN);
        P2SoftSphericalTruncated interpotentialNN = new P2SoftSphericalTruncated(space, new P2ElectrostaticDreiding(space, 2638.028500 * 11604.45728, 0.264550, 14.286224 * 11604.45728), truncationRadiusNN);
        potentialMaster.addPotential(interpotentialNN, new AtomType[]{species.getNType(), species.getNType()});

        ((CriterionInterMolecular) potentialMaster.getCriterion(interpotentialCC)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular) potentialMaster.getCriterion(interpotentialCO)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular) potentialMaster.getCriterion(interpotentialCN)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular) potentialMaster.getCriterion(interpotentialON)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular) potentialMaster.getCriterion(interpotentialOO)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular) potentialMaster.getCriterion(interpotentialNN)).setIntraMolecularCriterion(new CriterionNone());

        bdry = new BoundaryRectangularPeriodic(space, 1); //unit cell
        bdry.setBoxSize(space.makeVector(new double[]{2 * 17.248, 2 * 12.086, 4 * 7.382}));
        box.setBoundary(bdry);
        configOrthoLattice.initializeCoordinates(box);

        CoordinateDefinitionParacetamol coordDef = new CoordinateDefinitionParacetamol(this, box, primitive, basis, space);

    } //end of constructor

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        MDParacetamolOrthorhombic sim = new MDParacetamolOrthorhombic();
        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME,1);
        Pixel pixel = new Pixel(10);
        simGraphic.getDisplayBox(sim.box).setPixelUnit(pixel);
        ArrayList dataStreamPumps = simGraphic.getController().getDataStreamPumps();

        /*****************************************************************************/
        MeterKineticEnergy meterKE = new MeterKineticEnergy(sim.box);
        DisplayTextBox KEbox = new DisplayTextBox();
        DataPump KEpump = new DataPump(meterKE, KEbox);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(KEpump));
        dataStreamPumps.add(KEpump);

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster, sim.box);
        DisplayTextBox PEbox = new DisplayTextBox();
        DataPump PEpump = new DataPump(meterPE, PEbox);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(PEpump));
        dataStreamPumps.add(PEpump);

        MeterEnergy meterTotal = new MeterEnergy(sim.potentialMaster, sim.box);
        DisplayTextBox meterTotalbox = new DisplayTextBox();
        DataPump meterTotalpump = new DataPump(meterTotal, meterTotalbox);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterTotalpump));
        dataStreamPumps.add(meterTotalpump);

        MeterTemperature meterTemp = new MeterTemperature(sim.box, sim.space.D());
        DisplayTextBox tempBox = new DisplayTextBox();
        tempBox.setUnit(Kelvin.UNIT);
        DataPump tempPump = new DataPump(meterTemp, tempBox);
        sim.integrator.getEventManager().addListener(new IntegratorListenerAction(tempPump));
        dataStreamPumps.add(tempPump);


        /**********************************************************************/
        simGraphic.add(KEbox);
        simGraphic.add(PEbox);
        simGraphic.add(meterTotalbox);
        simGraphic.add(tempBox);


        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        colorScheme.setColor(sim.species.getAtomType(0), Color.red);
        colorScheme.setColor(sim.species.getAtomType(1), Color.gray);
        colorScheme.setColor(sim.species.getAtomType(2), Color.blue);
        colorScheme.setColor(sim.species.getAtomType(3), Color.white);
        colorScheme.setColor(sim.species.getAtomType(4), Color.white);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        simGraphic.getDisplayBox(sim.box).repaint();

    }//end of main
}//end of class
