package etomica.models.clathrates;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.SimulationGraphic;
import etomica.models.clathrates.MinimizationTIP4P.ChargeAgentSourceRPM;
import etomica.models.water.ConfigurationFileTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMoleculeList;
import etomica.potential.EwaldSummation;
import etomica.potential.EwaldSummation.MyCharge;
import etomica.potential.P2LennardJones;
import etomica.potential.Potential2SoftSpherical;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Calorie;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.util.Constants;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Implementation of Buch's algorithm for proton disorder.
 * <p>
 * Buch, V.; Sandler, P.; Sadlej, J. Simulations of H2O Solid, Liquid,
 * and Clusters, with an Emphasis on Ferroelectric Ordering Transition
 * in Hexagonal Ice. J. Phys. Chem. B 1998, 102, 8641−8653.
 *
 * @author Sabry Moustafa
 */
public class ProtonDisorderGenerator extends Simulation {
    protected Box boxO, box;
    protected SpeciesSpheresMono speciesO;
    protected SpeciesWater4P species;
    protected Potential2SoftSpherical potentialLJ;
    protected EwaldSummation potentialES;
    protected PotentialMaster potentialMaster;

    public ProtonDisorderGenerator(Space space, String configFile, int nBasis, int[] nC, double[] a0, boolean isIce, int numMolecule, double rCutRealES, double kCut) {
        super(space);
        speciesO = new SpeciesSpheresMono(this, space);
        addSpecies(speciesO);
        Boundary boundaryO = new BoundaryRectangularPeriodic(space, a0);
        boxO = new Box(boundaryO, space);
        addBox(boxO);
        boxO.setNMolecules(speciesO, nBasis);

        ConfigurationFile config = new ConfigurationFile(configFile);
        config.initializeCoordinates(boxO);

        species = new SpeciesWater4P(space);
        addSpecies(species);
        double[] a0_sc = new double[]{a0[0] * nC[0], a0[1] * nC[1], a0[2] * nC[2]};
        Boundary boundary = new BoundaryRectangularPeriodic(space, a0_sc);
        box = new Box(boundary, space);
        addBox(box);
        box.setNMolecules(species, numMolecule);

//		double precision = 1.0e-5 , precision_s;
//		if (precision ==1.0e-5 ){   
//			precision_s = 3.047059472445871 ;
//		}	else if (precision == 5.0e-5){
//			precision_s = 2.800672811371045;}
//		else {throw new RuntimeException("improper precision value!");}

        ChargeAgentSourceRPM agentSource = new MinimizationTIP4P.ChargeAgentSourceRPM(species, isIce);
//		AtomLeafAgentManager atomAgentManager = new AtomLeafAgentManager(agentSource, box);
        AtomLeafAgentManager<MyCharge> atomAgentManager = new AtomLeafAgentManager<>(agentSource, box);

        double sigma, epsilon;
//		if(isIce){
//			sigma = 3.1668; epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice			
//		}else{			
//			sigma = 3.154; epsilon = Kelvin.UNIT.toSim(78.0); //TIP4P
//		}

        if (isIce) {
            sigma = 3.1668;
            epsilon = Kelvin.UNIT.toSim(106.1);//TIP4P/Ice
        } else {//TIP4P
            double A = 600E3; // kcal A^12 / mol
            double C = 610.0; // kcal A^6 / mol
            double s6 = A / C;
            sigma = Math.pow(s6, 1.0 / 6.0);
            epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000)) / 4.0;
        }

        potentialLJ = new P2LennardJones(space, sigma, epsilon);
        potentialLJ.setBox(box);
        potentialES = new EwaldSummation(box, atomAgentManager, space, kCut, rCutRealES);

//XXXX Potential Master
        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(potentialLJ, new AtomType[]{species.getOxygenType(), species.getOxygenType()});
        potentialMaster.addPotential(potentialES, new AtomType[0]);
    }

    public static void main(String[] args) {
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        String configFileName = params.configFile;
        int[] nC = params.nC;
        int nBasis = params.nBasis;
        double[] a0 = params.a0;

        double rCutRealES = params.rCutRealES;
        double kCut = params.kCut;

        boolean isIce = params.isIce;
        int numMolecule = nBasis * nC[0] * nC[1] * nC[2];
        int sample = params.sample;
        boolean includeM = params.includeM;


        boolean isGraphics = params.isGraphics;
        final ProtonDisorderGenerator sim = new ProtonDisorderGenerator(Space3D.getInstance(3), configFileName, nBasis, nC, a0, isIce, numMolecule, rCutRealES, kCut);
        int[][] pairOdCoord = new int[nBasis * 4 / 2][3];
        int[] nCoordsO = new int[nBasis];
        IAtomList atomList = sim.boxO.getLeafList();
        Vector dr = sim.space.makeVector();
        double dr2, drNbr = 3.5;
        int nCoordTot = 0; //must = N*4/2
        for (int i = 0; i < nBasis; i++) {
            for (int j = i + 1; j < nBasis; j++) {
                if (i == j) continue;
                dr.Ev1Mv2(atomList.get(j).getPosition(), atomList.get(i).getPosition());
                sim.boxO.getBoundary().nearestImage(dr);
                dr2 = dr.squared();
                if (dr2 < drNbr * drNbr) {
                    int dRandCoordi = (sim.getRandom().nextDouble() < 0.5) ? 0 : 1;
                    int dRandCoordj = 1 - dRandCoordi;
                    nCoordsO[i] += dRandCoordi;
                    nCoordsO[j] += dRandCoordj;
                    pairOdCoord[nCoordTot][0] = i;
                    pairOdCoord[nCoordTot][1] = j;
                    pairOdCoord[nCoordTot][2] = dRandCoordi;
                    nCoordTot++;
                }
            }//j
        }//i
        int test = 0, count = 0;
        while (test != nBasis) {
            int H = sim.getRandom().nextInt(nCoordTot);
            int iO = pairOdCoord[H][0];
            int jO = pairOdCoord[H][1];
            int dCoordNew;

            int dCoordOld = Math.abs(nCoordsO[iO] - nCoordsO[jO]);
            if (pairOdCoord[H][2] == 1) {
                dCoordNew = Math.abs((nCoordsO[iO] - 1) - (nCoordsO[jO] + 1));
            } else {
                dCoordNew = Math.abs((nCoordsO[iO] + 1) - (nCoordsO[jO] - 1));
            }
            if (dCoordNew < dCoordOld) {
                nCoordsO[iO] += (pairOdCoord[H][2] == 1 ? -1 : 1);
                nCoordsO[jO] += (pairOdCoord[H][2] == 1 ? 1 : -1);
                pairOdCoord[H][2] = (pairOdCoord[H][2] == 1) ? 0 : 1;
            } else if (dCoordNew == dCoordOld && sim.getRandom().nextDouble() < 0.5) {
                nCoordsO[iO] += (pairOdCoord[H][2] == 1 ? -1 : 1);
                nCoordsO[jO] += (pairOdCoord[H][2] == 1 ? 1 : -1);
                pairOdCoord[H][2] = (pairOdCoord[H][2] == 1) ? 0 : 1;
            }//End if
            test = 0;
            count++;
            for (int k = 0; k < nBasis; k++) {
                if (nCoordsO[k] == 2) {
                    test++;
                }
            }
            if (test == nBasis) {
                System.out.println("Converged in " + count + " MC steps");
                break;
            }
        }//End MC

        int[][] OH12 = new int[nBasis][2];
        for (int i = 0; i < nBasis; i++) {
            OH12[i][0] = -1;
            OH12[i][1] = -1;
        }
        for (int i = 0; i < pairOdCoord.length; i++) {//h
            int Oi = pairOdCoord[i][1 - pairOdCoord[i][2]];
            if (OH12[Oi][0] == -1) {
                OH12[Oi][0] = i;
            } else if (OH12[Oi][1] == -1) {
                OH12[Oi][1] = i;
            } else {
                throw new RuntimeException();
            }
        }
        double rOH = 0.9572;
        Vector posH = sim.space.makeVector();
        IMoleculeList molList = sim.box.getMoleculeList();
        String confProtonDisorder = "config_" + sample;
        try {
            PrintWriter outputStreamM = new PrintWriter(confProtonDisorder + ".pos");
            for (int i = 0; i < nBasis; i++) {
                Vector posO0 = atomList.get(i).getPosition();
                molList.get(i).getChildList().get(0).getPosition().E(posO0);
                for (int H = 0; H < 2; H++) {
                    int j = (pairOdCoord[OH12[i][H]][0] == i) ? pairOdCoord[OH12[i][H]][1] : pairOdCoord[OH12[i][H]][0];
                    Vector posO1 = atomList.get(j).getPosition();
                    dr.Ev1Mv2(posO1, posO0);
                    sim.boxO.getBoundary().nearestImage(dr);
                    dr.normalize();
                    posH.E(posO0);
                    posH.PEa1Tv1(rOH, dr);
                    molList.get(i).getChildList().get(H + 1).getPosition().E(posH);
                    outputStreamM.println(posH.getX(0) + " " + posH.getX(1) + " " + posH.getX(2));
                }
                if (i == (nBasis - 1)) {
                    outputStreamM.print(atomList.get(i).getPosition().getX(0) + " " + atomList.get(i).getPosition().getX(1) + " " + atomList.get(i).getPosition().getX(2));
                } else {
                    outputStreamM.println(atomList.get(i).getPosition().getX(0) + " " + atomList.get(i).getPosition().getX(1) + " " + atomList.get(i).getPosition().getX(2));
                }
            }//i->nBasis
            outputStreamM.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }

        if (includeM) {
            ConfigurationFile config = new ConfigurationFile(configFileName);////to duplicate with M point!
            ConfigurationFileBinary.replicate(config, sim.box, nC, sim.space);
        } else {
            ConfigurationFileTIP4P config = new ConfigurationFileTIP4P(confProtonDisorder, sim.space, isIce);
            ConfigurationFileBinary.replicate(config, sim.box, nC, sim.space);
        }


//        config.initializeCoordinates(sim.box);
        Vector dipoleP = sim.space.makeVector();
        for (int i = 0; i < numMolecule; i++) {
            for (int j = 0; j < 2; j++) {
                dr.Ev1Mv2(molList.get(i).getChildList().get(j).getPosition(), molList.get(i).getChildList().get(3).getPosition());
                sim.boxO.getBoundary().nearestImage(dr);
                dipoleP.PEa1Tv1(0.52, dr);
            }
//		    System.out.println(Math.sqrt(dipoleP.squared())/0.20819434); // units in D
        }
        System.out.println("P = " + Math.sqrt(dipoleP.squared()) / 0.20819434 / 2.18 + " (/2.18D) "); // in "2.18 D" units


        MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(sim.potentialMaster);
        meterPotentialEnergy.setBox(sim.box);
        double E = Joule.UNIT.fromSim(meterPotentialEnergy.getDataAsScalar() / sim.box.getMoleculeList().size()) * 1.0E-3 * Constants.AVOGADRO;
        System.out.println("E (kJ/mol)  = " + E);

        if (isGraphics) {
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, "string");
            final DisplayBox display = new DisplayBox(sim, sim.box);
            simGraphic.add(display);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getHydrogenType(), Color.GREEN);
            ((ColorSchemeByType) simGraphic.getDisplayBox(sim.box).getColorScheme()).setColor(sim.species.getOxygenType(), Color.RED);
            //Sabry
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getOxygenType(), 2.0);
            ((DiameterHashByType) ((DisplayBox) simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getHydrogenType(), 1.0);

            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBackgroundColor(Color.white);
            ((DisplayBoxCanvasG3DSys) simGraphic.getDisplayBox(sim.box).canvas).setBoundaryFrameColor(Color.blue);
            simGraphic.makeAndDisplayFrame();
        }


    }//END main

    public static class SimParams extends ParameterBase {
        public int nBasis = 68;

        public String configFile = "configOxygenSH";
        int nX = 1;
        public int[] nC = new int[]{nX, nX, nX};
        //		public double[] a0 = new double[]{12.03, 12.03, 12.03};//sI
        public double[] a0 = new double[]{12.21, 21.15, 10.14};//sI
        public double rCutRealES = 20.0;
        public double kCut = 2.0;
        public boolean isIce = false, isGraphics = false;
        public int sample = 1;
        public boolean includeM = false;

    }
}