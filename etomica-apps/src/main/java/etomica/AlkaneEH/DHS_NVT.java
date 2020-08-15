/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.AlkaneEH;

import etomica.action.BoxImposePbc;

import etomica.action.activity.ActivityIntegrate2;
import etomica.action.controller.Controller;
import etomica.atom.AtomTypeOriented;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPump;
import etomica.data.meter.MeterDipoleSumSquared1site;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayBoxCanvasG3DSys.OrientedFullSite;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotate;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.potential.P2HSDipole;
import etomica.potential.P2MoleculeTruncated;
import etomica.potential.P2ReactionFieldDipole;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;

/**
 * Canonical ensemble Monte Carlo simulation (NVT)
 * calculate dielectric constant for dipolar HS model
 * 
 * @author shu
 * Dec 2014
 */

public class DHS_NVT extends Simulation {
	private static final long serialVersionUID = 1L;
    private final static String APP_NAME = "dipolar HS, dielectric constant";
    private static final int PIXEL_SIZE = 15;
    protected final SpeciesSpheresRotating species;
    protected final PotentialMaster potentialMaster;
	protected final IntegratorMC integrator;
	protected final MCMoveMolecule moveMolecule;//translation mc move
	protected final MCMoveRotate rotateMolecule;//atomic rotation mc move
	protected final Box box;

	//************************************* constructor ********************************************//
    public DHS_NVT(Space space, int numberMolecules, final double sigmaHS, double mu,
                   double dielectricOutside, double boxSize, double temperature, double truncation) {
        super(space);
        //setRandom(new RandomNumberGenerator(1));
        species = new SpeciesSpheresRotating(space, new ElementSimple("A"));
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numberMolecules);
        box.getBoundary().setBoxSize(Vector.of(new double[]{boxSize, boxSize, boxSize}));

        IMoleculePositionDefinition positionDefinition = new IMoleculePositionDefinition() {
            public Vector position(IMolecule molecule) {
                return molecule.getChildList().get(0).getPosition();
            }
        };

        // potential part
        P2HSDipole pTarget = new P2HSDipole(space, sigmaHS, mu);
        DipoleSourceDHS dipoleDHS = new DipoleSourceDHS(space, mu);// add reaction field potential
        System.out.println("in main class, magnitude of dipole:" + dipoleDHS.dipoleStrength);
        P2ReactionFieldDipole pRF = new P2ReactionFieldDipole(space, positionDefinition);
        pRF.setDipoleSource(dipoleDHS);
        pRF.setRange(truncation);
        pRF.setDielectric(dielectricOutside);

        potentialMaster = new PotentialMaster();
        P2MoleculeTruncated p2TruncatedRF = new P2MoleculeTruncated(pRF, truncation, space, positionDefinition);
        potentialMaster.addPotential(p2TruncatedRF, new ISpecies[]{species, species});//add truncated potential from reaction field to potential master
        potentialMaster.lrcMaster().addPotential(pRF.makeP0());// add P0ReactionField potential to potential master
        // u(HS)+u(dd)
        P2MoleculeTruncated p2TruncatedDHS = new P2MoleculeTruncated(pTarget, truncation, space, positionDefinition);
        potentialMaster.addPotential(p2TruncatedDHS, new ISpecies[]{species, species});

        integrator = new IntegratorMC(this, potentialMaster, box);
        moveMolecule = new MCMoveMolecule(this, potentialMaster, space);
        rotateMolecule = new MCMoveRotate(potentialMaster, random, space);

        this.getController2().addActivity(new ActivityIntegrate2(integrator));

        //******************************** periodic boundary condition ******************************** //
        BoxImposePbc imposePbc = new BoxImposePbc(box, space);
        imposePbc.setApplyToMolecules(true);

        integrator.setTemperature(temperature);
        integrator.getMoveManager().addMCMove(moveMolecule);
        integrator.getMoveManager().addMCMove(rotateMolecule);
        integrator.getEventManager().addListener(new IntegratorListenerAction(imposePbc));

        //******************************** initial configuration ******************************** //
        LatticeCubicFcc lattice = new LatticeCubicFcc(space);
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
		long t1 = System.currentTimeMillis();
		Space space = Space3D.getInstance();
		int steps = params.steps;
		boolean isGraphic = params.isGraphic;
		double temperature = params.temperature;//in sim unit
		int numberMolecules = params.numberMolecules;
		double density = params.density;
		double sigmaHS=params.sigmaHS;
		double dipoleStrength = Math.sqrt(params.dipoleStrength2);
		double dielectricOutside = params.dielectricOutside;
		double densitySim = density;
		double boxSize = Math.pow(numberMolecules/densitySim,(1.0/3.0)); //unit of angstrom
		double truncation=boxSize * 0.49;

		System.out.println("******************* dipolar HS, dielectric constant, NVT********************");
		System.out.println("steps= "+steps);
		System.out.println("number of molecules= "+numberMolecules);
		System.out.println("temperature(sim) ="+temperature);
		System.out.println("density(sim)1/angstrom^3 ="+densitySim);
		System.out.println("box size(angstrom):"+boxSize);
		System.out.println("truncation="+truncation);
		System.out.println("sigmaHS="+sigmaHS);
		System.out.println("dipoleStrength="+dipoleStrength);
		System.out.println("dielectricOutside="+dielectricOutside);

		final DHS_NVT sim = new DHS_NVT(space,numberMolecules,sigmaHS,dipoleStrength,dielectricOutside,boxSize,temperature,truncation);

    	if (isGraphic){

			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
		    simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(PIXEL_SIZE));
            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
            ((DiameterHashByType)((DisplayBox)simGraphic.displayList().getFirst()).getDiameterHash()).setDiameter(sim.species.getAtomType(0),1);
            ColorSchemeByType colorScheme = (ColorSchemeByType)simGraphic.getDisplayBox(sim.box).getColorScheme();
            colorScheme.setColor(sim.getSpecies(0).getAtomType(0), Color.red);
            OrientedFullSite[] sites = new OrientedFullSite[2];
            sites[0] = new OrientedFullSite(Vector.of(new double[]{0.5, 0, 0}), Color.BLUE, 0.2);
            sites[1] = new OrientedFullSite(Vector.of(new double[]{-0.5, 0, 0}), Color.YELLOW, 0.2);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setOrientationSites(
                    (AtomTypeOriented) sim.getSpecies(0).getAtomType(0), sites);
            ((DisplayBoxCanvasG3DSys)simGraphic.getDisplayBox(sim.box).canvas).setOrientationSites(
                    (AtomTypeOriented) sim.getSpecies(0).getAtomType(0), sites);
            simGraphic.makeAndDisplayFrame(APP_NAME);
			simGraphic.getDisplayBox(sim.box).repaint();
	    	return ;
    	}
    	sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), steps/5);// equilibration period

   		sim.integrator.getMoveManager().setEquilibrating(false);
   		System.out.println("equilibration finished");

        int blockNumber = 100;
   		int sampleAtInterval = numberMolecules;
   		int samplePerBlock = steps/sampleAtInterval/blockNumber;
   		System.out.println("number of blocks is : "+blockNumber);
   		System.out.println("sample per block is : "+samplePerBlock);

        // dipoleSumSquared
        MeterDipoleSumSquared1site dipoleSumSquaredMeter = new MeterDipoleSumSquared1site(space, sim.box, dipoleStrength);
        AccumulatorAverage dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
        DataPump dipolePump = new DataPump(dipoleSumSquaredMeter,dipoleSumSquaredAccumulator);
        IntegratorListenerAction dipoleListener = new IntegratorListenerAction(dipolePump);
        dipoleListener.setInterval(sampleAtInterval);
        sim.integrator.getEventManager().addListener(dipoleListener);
        // energy
        MeterPotentialEnergyFromIntegrator energyMeter = new MeterPotentialEnergyFromIntegrator(sim.integrator);
        AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
        DataPump energyPump = new DataPump(energyMeter, energyAccumulator);
        energyAccumulator.setBlockSize(50);
        IntegratorListenerAction energyListener = new IntegratorListenerAction(energyPump);
        sim.integrator.getEventManager().addListener(energyListener);
sim.getController2().runActivityBlocking(new ActivityIntegrate2(sim.integrator), steps);

        //calculate dipoleSumSquared average
        double dipoleSumSquared = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.AVERAGE.index)).x;
        double dipoleSumSquaredERR = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.ERROR.index)).x;
        double volume = sim.box.getBoundary().volume();
        double dipoleFac = 4 * Math.PI * dipoleSumSquared/9.0/volume/temperature;
        double dielectricOutsideFac = 2*(dielectricOutside-1)/(2*dielectricOutside+1);
        double x1 =  dipoleSumSquared;
        double B  = dielectricOutsideFac;
        double D  = 4 * Math.PI/9.0/volume/temperature;
        System.out.println("B:  "+B);
        System.out.println("D:  "+D);
        System.out.println("x1:  "+x1);

        double dEpsilondx_1=(B*D+2*D)/(B*D*x1-D*x1+1);
        double dEpsilondx_2=(B*D-D)*(B*D*x1+2D*x1+1)/Math.pow((B*D*x1-D*x1+1),2);
        double dEpsilondx =dEpsilondx_1-dEpsilondx_2;
        double epsilon = (1+B*D*x1+2*D*x1)/(1+B*D*x1-D*x1);
        double epsilonERR_alt = dEpsilondx*dipoleSumSquaredERR;
        double epsilonERR = 3 * D / Math.pow((B*D*x1-D*x1+1),2)*dipoleSumSquaredERR;

        System.out.println("epsilonERR_alt: "+epsilonERR_alt);
        System.out.println("<M^2> :  "+dipoleSumSquared + "   with err:  "+dipoleSumSquaredERR );
        System.out.println("Epsilon="+epsilon + "   with err: "+epsilonERR);
        System.out.println("========================");
        double C = dipoleFac;
        double A = C/(1+B*C);
        double dielectricConstant = (1+2*A)/(1-A);
        System.out.println("(epsilon-1)/(epsilon+2): "+A);
        System.out.println("dielectric constant is:  "+dielectricConstant);

        double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index)).x;
        avgPE /= numberMolecules;
        System.out.println("PE/epsilon:"+avgPE);
        long t2 = System.currentTimeMillis();
        System.out.println("simulation time is:"+ (t2-t1));
	}

    public static class DipoleSourceDHS implements DipoleSource {//for potential reaction field
        protected final Vector dipoleVector;
        protected double dipoleStrength;

        public DipoleSourceDHS(Space space, double dipole) {
            dipoleStrength = dipole;
            dipoleVector = space.makeVector();
        }

        public Vector getDipole(IMolecule molecule) {
            if (molecule.getChildList().size() != 1) {
                throw new RuntimeException("improper number of atom in the molecule");
            }
            IAtomList atomList = molecule.getChildList();
            IAtomOriented atom = (IAtomOriented) atomList.get(0);
            dipoleVector.E(atom.getOrientation().getDirection());
            dipoleVector.TE(dipoleStrength);
            return dipoleVector;
        }
    }

    // ******************* parameters **********************//
	public static class Param extends ParameterBase {
		public boolean isGraphic = false;
		public double temperature = 2.0;
		public int numberMolecules = 100;
		public double density =0.5; 
		public double dipoleStrength2=1.0;
		public double sigmaHS=1.0;
		public double dielectricOutside = 1000;
		public int steps = 100000;
	}
}
