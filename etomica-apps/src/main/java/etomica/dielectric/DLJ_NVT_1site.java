/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dielectric;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeOriented;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterDipoleSumSquared1site;
import etomica.data.meter.MeterDipoleSumSquaredMappedAverage;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayBoxCanvasG3DSys.OrientedSite;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicBcc;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.potential.P2LJDipole;
import etomica.potential.P2ReactionFieldDipole;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Canonical ensemble Monte Carlo simulation (NVT)
 * calculate dielectric constant epsilon for dipolar LJ model
 * @author Weisong & shu
 * July 2015
 */

public class DLJ_NVT_1site extends Simulation {

	private static final long serialVersionUID = 1L;
    private final static String APP_NAME = "dipolar LJ";
    private static final int PIXEL_SIZE = 15;
    protected final PotentialMaster potentialMaster;
	protected final IntegratorMC integrator;
	protected final MCMoveMolecule moveMolecule;//translation
	protected final MCMoveRotate rotateMolecule;//rotation, atomic
	protected final Box box;
    protected SpeciesGeneral species;

	//************************************* constructor ********************************************//
    public DLJ_NVT_1site(Space space, int numberMolecules, final double sigmaLJ, double epsilonLJ, double mu,
                         double dielectricOutside, double boxSize, double temperature, double truncation) {
        super(space);
        species = SpeciesSpheresRotating.create(space, new ElementSimple("A"));
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numberMolecules);
        box.getBoundary().setBoxSize(Vector.of(new double[]{boxSize, boxSize, boxSize}));

        IMoleculePositionDefinition positionDefinition = new IMoleculePositionDefinition() {
            public Vector position(IMolecule molecule) {
                return molecule.getChildList().get(0).getPosition();
            }
        };

        // dipolar LJ potential
        P2LJDipole pDLJ = new P2LJDipole(space, sigmaLJ, epsilonLJ, mu, truncation);
        // add reaction field potential
        DipoleSourceDLJ dipoleSourceDLJ = new DipoleSourceDLJ(space, mu);// add reaction field potential
        P2ReactionFieldDipole pRF = new P2ReactionFieldDipole(space, positionDefinition);
        pRF.setDipoleSource(dipoleSourceDLJ);
        pRF.setRange(truncation);
        pRF.setDielectric(dielectricOutside);


        potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(pDLJ, new ISpecies[]{species, species});
        potentialMaster.addPotential(pRF, new ISpecies[]{species, species});
        potentialMaster.lrcMaster().addPotential(pRF.makeP0());

        // integrator from potential master
        integrator = new IntegratorMC(this, potentialMaster, box);
        // add mc move
        moveMolecule = new MCMoveMolecule(this, potentialMaster, space);//stepSize:1.0, stepSizeMax:15.0  ??????????????
        rotateMolecule = new MCMoveRotate(potentialMaster, random, space);

        this.getController().addActivity(new ActivityIntegrate(integrator));

        //******************************** periodic boundary condition ******************************** //
        BoxImposePbc imposePbc = new BoxImposePbc(box, space);
        imposePbc.setApplyToMolecules(true);
        //**************************** integrator ****************************** //
        integrator.setTemperature(temperature);
        integrator.getMoveManager().addMCMove(moveMolecule);
        integrator.getMoveManager().addMCMove(rotateMolecule);
        integrator.getEventManager().addListener(new IntegratorListenerAction(imposePbc));

        //******************************** initial configuration ******************************** //
        LatticeCubicBcc lattice = new LatticeCubicBcc(space);
        ConfigurationLattice configuration = new ConfigurationLattice(lattice, space);
        configuration.initializeCoordinates(box);

    }

    // **************************** simulation part **************************** //
	public static void main (String[] args){
		Param params = new Param();
        if (args.length > 0) {
            ParseArgs.doParseArgs(params, args);
		} else {

		}
		final long startTime = System.currentTimeMillis();
		Space space = Space3D.getInstance();
		int steps = params.steps;
		boolean isGraphic = params.isGraphic;
		boolean mSquare = params.mSquare;
		boolean aEE = params.aEE;

		double temperature = params.temperature;
	//	double temperatureK = params.temperatureK;
//		double temperature = Kelvin.UNIT.toSim(temperatureK);// convert Kelvin temperature to T(sim), essentially kT
		int numberMolecules = params.numberMolecules;
		double density = params.density;
		double sigmaLJ=params.sigmaLJ;
		double epsilonLJ=params.epsilonLJ;
		double dipoleStrength = Math.sqrt(params.dipoleStrength2);
		double dielectricOutside = params.dielectricOutside;
		//double densitySim = density * Constants.AVOGADRO * 1e-27;  // ?????? convert density to sim unit; in 1/(A)^3
		//System.out.println("Constants.AVOGADRO * 1e-27: "+Constants.AVOGADRO * 1e-27);
		double boxSize = Math.pow(numberMolecules/density,(1.0/3.0));
		double truncation=boxSize* 0.49;


//		System.out.println("******************* dipolar LJ, dielectric constant, NVT********************");
		System.out.println("number of molecules =\t"+numberMolecules);
		System.out.println("steps= "+ steps);
		System.out.println("density=\t"+density);
//		System.out.println("denisty(sim)="+densitySim);
		System.out.println("temperature=\t"+ temperature);
//		System.out.println("box size="+boxSize);
//		System.out.println("truncation = "+truncation);
//		System.out.println("sigmaLJ="+sigmaLJ);
//		System.out.println("epsilonLJ="+epsilonLJ);
		System.out.println("dipoleStrength squared = "+params.dipoleStrength2);
//		System.out.println("dipoleStrength="+dipoleStrength);
		System.out.println("dielectricOutside="+dielectricOutside);

		final DLJ_NVT_1site sim = new DLJ_NVT_1site(space,numberMolecules,sigmaLJ, epsilonLJ, dipoleStrength,
				dielectricOutside, boxSize,temperature,truncation);

    	if (isGraphic){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getAtomType(0),1);
            ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.red);
            OrientedSite[] sites = new OrientedSite[1];
            sites[0] = new OrientedSite(0.5, Color.BLUE, 0.2);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setOrientationSites(
                    (AtomTypeOriented) sim.getSpecies(0).getAtomType(0), sites);
            simGraphic.makeAndDisplayFrame(APP_NAME);
			simGraphic.getDisplayBox(sim.box).repaint();
	    	return ;
    	}

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps / 5));// equilibration period

   		sim.integrator.getMoveManager().setEquilibrating(false);
//   		System.out.println("equilibration finished");

        int blockNumber = 1000;
   		int sampleAtInterval = numberMolecules;
   		int samplePerBlock = steps/sampleAtInterval/blockNumber;
//   		System.out.println("number of blocks is : "+blockNumber);
//   		System.out.println("sample per block is : "+samplePerBlock);

        MeterDipoleSumSquared1site dipoleSumSquaredMeter = null;
   		AccumulatorAverage dipoleSumSquaredAccumulator = null;
   		if(mSquare){
            dipoleSumSquaredMeter = new MeterDipoleSumSquared1site(space, sim.box, dipoleStrength);
            dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
   			DataPump dipolePump = new DataPump(dipoleSumSquaredMeter,dipoleSumSquaredAccumulator);
   			IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
   			dipoleListener.setInterval(sampleAtInterval);
   			sim.integrator.getEventManager().addListener(dipoleListener);
   		}

        //AEE
        DipoleSourceDLJ dipoleSourceDLJ = new DipoleSourceDLJ(space, dipoleStrength);
        MeterDipoleSumSquaredMappedAverage AEEMeter =  null;
   		AccumulatorAverageCovariance AEEAccumulator = null;
   		if(aEE){
         AEEMeter = new MeterDipoleSumSquaredMappedAverage(sim.box, sim.getSpeciesManager(), dipoleStrength, temperature, sim.potentialMaster);
            AEEMeter.setDipoleSource(dipoleSourceDLJ);
         AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock,true);
		DataPump AEEPump = new DataPump(AEEMeter,AEEAccumulator);
		IntegratorListenerAction AEEListener = new IntegratorListenerAction(AEEPump);


            AEEListener.setInterval(sampleAtInterval);
//		AEEListener.setInterval(1);// debug only

            sim.integrator.getEventManager().addListener(AEEListener);

        }
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

        //calculate dipoleSumSquared average
        double dipoleSumSquared = 0;
        double dipoleSumSquaredERR = 0;
        double dipoleSumCor = 0 ;
        if(mSquare){
            dipoleSumSquared = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.AVERAGE.index)).x;
            dipoleSumSquaredERR = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.ERROR.index)).x;
            dipoleSumCor = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.BLOCK_CORRELATION.index)).x;
        }

        double AEE = 0;
        double AEEER =0;
        double AEECor = 0;
        if(aEE){
            double sum0 = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.AVERAGE.index).getValue(0);
            double ERsum0 = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(0);
            AEECor = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_CORRELATION.index).getValue(0);
            AEE = sum0;
		 AEEER = ERsum0;
        }

        long endTime = System.currentTimeMillis();

        double totalTime = (endTime - startTime)/1000.0;
        if(mSquare){
        System.out.println("-<M^2>*bt*bt:\t"+(-dipoleSumSquared/temperature/temperature)
        		+ " mSquareErr:\t" + (dipoleSumSquaredERR/temperature/temperature)
        		+ " mSquareDifficulty:\t"+(dipoleSumSquaredERR/temperature/temperature)*Math.sqrt(totalTime)
        		+ " dipolesumcor= " + dipoleSumCor );
        }
        if(aEE){
            System.out.println("AEE_new:\t" + (AEE)
                    + " AEEErr:\t" + AEEER
                    + " AEEDifficulty:\t"+ AEEER*Math.sqrt(totalTime)
        		+ " AEECor= " + AEECor );
        }

        System.out.println(  "time: " + totalTime);
	}

    public static class DipoleSourceDLJ implements DipoleSource {//for potential reaction field
        protected final Vector dipoleVector;
        protected double dipoleStrength;

        public DipoleSourceDLJ(Space space, double s) {
            dipoleStrength = s;
            dipoleVector = space.makeVector();
        }

        public Vector getDipole(IMolecule molecule) {
            IAtomList atomList = molecule.getChildList();
            if (atomList.size() != 1) {
                throw new RuntimeException("improper input of the molecule");
            }
            IAtomOriented atom = (IAtomOriented) atomList.get(0);
            dipoleVector.E(atom.getOrientation().getDirection());
            dipoleVector.TE(dipoleStrength);
            return dipoleVector;
        }
    }

    // ******************* parameters **********************//
	public static class Param extends ParameterBase {
		public boolean isGraphic = false;
		public boolean mSquare = false;
		public boolean aEE = true;
		public double temperature = 10;
		public int numberMolecules = 10;
		public double density = 0.325;
		public double dipoleStrength2 = 3;
		public double sigmaLJ = 1.0;
		public double epsilonLJ = 1.0;
		public double dielectricOutside = 1.0E11;
		public int steps = 1000000;
	}
}
