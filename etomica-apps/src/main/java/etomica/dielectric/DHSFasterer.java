/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dielectric;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomTypeOriented;
import etomica.atom.DiameterHashByType;
import etomica.atom.IAtom;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.data.*;
import etomica.data.meter.MeterDipoleSumSquared1site;
import etomica.data.meter.MeterDipoleSumSquaredMappedAverageFasterer;
import etomica.data.meter.MeterPotentialEnergyFromIntegratorFasterer;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxCanvasG3DSys;
import etomica.graphics.DisplayBoxCanvasG3DSys.OrientedSite;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveAtomFasterer;
import etomica.integrator.mcmove.MCMoveAtomRotateFasterer;
import etomica.lattice.LatticeCubicFcc;
import etomica.molecule.DipoleSourceAtomic;
import etomica.potential.P1ReactionField;
import etomica.potential.P2HSDipoleAtomic;
import etomica.potential.P2ReactionFieldDipoleFasterer;
import etomica.potential.P2SoftSum;
import etomica.potential.compute.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.SpeciesGeneral;
import etomica.species.SpeciesSpheresRotating;
import etomica.units.Pixel;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

/**
 * Canonical ensemble Monte Carlo simulation (NVT)
 * calculate dielectric constant for dipolar HS model
 * 
 * Electric field mapped average for dielectric constant
 * 
 * @author Shu & Weisong
 * May 2015
 */

public class DHSFasterer extends Simulation {
    private final static String APP_NAME = "dipolar HS, dielectric constant";
    private static final int PIXEL_SIZE = 15;
	protected final SpeciesGeneral species;
    protected final PotentialCompute potentialMaster;
	protected final IntegratorMCFasterer integrator;
	protected final MCMoveAtomFasterer moveMolecule;//translation mc move
	protected final MCMoveAtomRotateFasterer rotateMolecule;//atomic rotation mc move
	protected final Box box;

	//************************************* constructor ********************************************//
    public DHSFasterer(Space space, int numberMolecules, final double HSDiameter, double mu,
                       double dielectricOutside, double boxSize, double temperature, double truncation) {
        super(space);

        species = SpeciesSpheresRotating.create(space, new ElementSimple("A"), false,true);
        addSpecies(species);
        box = this.makeBox();
        box.setNMolecules(species, numberMolecules);
        box.getBoundary().setBoxSize(Vector.of(new double[]{boxSize, boxSize, boxSize}));

        P2HSDipoleAtomic pTarget = new P2HSDipoleAtomic(space, HSDiameter, mu, truncation);
        DipoleSourceDHS dipoleDHS = new DipoleSourceDHS(space, mu);// add reaction field potential
//		System.out.println("in main class, magnitude of dipole:"+dipoleDHS.dipoleStrength);
        P2ReactionFieldDipoleFasterer pRF = new P2ReactionFieldDipoleFasterer(space, dipoleDHS);
        pRF.setRange(truncation);
        pRF.setDielectric(dielectricOutside);

        P2SoftSum p2 = new P2SoftSum(pTarget, pRF);

		NeighborManagerSimple neighborManager = new NeighborManagerSimple(box);
        PotentialComputePairGeneral potentialPair = new PotentialComputePairGeneral(getSpeciesManager(), box, neighborManager);

        potentialPair.setPairPotential(species.getLeafType(), species.getLeafType(), p2);

		PotentialComputeField pcField = new PotentialComputeField(getSpeciesManager(), box);
		P1ReactionField p1RF = new P1ReactionField(dipoleDHS, dielectricOutside, truncation);
		pcField.setFieldPotential(species.getLeafType(), p1RF);

		potentialMaster = new PotentialComputeAggregate(pcField, potentialPair);

        // integrator from potential master
        integrator = new IntegratorMCFasterer(potentialMaster, random, temperature, box);
        moveMolecule = new MCMoveAtomFasterer(random, potentialMaster, box);
		rotateMolecule = new MCMoveAtomRotateFasterer(random, potentialMaster, box);

        this.getController().addActivity(new ActivityIntegrate(integrator));

        //******************************** periodic boundary condition ******************************** //
        BoxImposePbc imposePbc = new BoxImposePbc(box, space);
        imposePbc.setApplyToMolecules(true);
        //**************************** integrator ****************************** //
        integrator.setTemperature(temperature);
        integrator.getMoveManager().addMCMove(moveMolecule);  //TODO
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
		    params.density = 0.0001;
		    params.dielectricOutside = 1.00209587222; // 0.001
		    params.dielectricOutside = 1.00020945428; // 0.0001
		    params.dipoleStrength2 = 0.5;
		    params.temperature = 1;
		    params.numberMolecules = 2;
		    params.steps = 1000000000L;
		    params.isGraphic = false;
		}
		final long startTime = System.currentTimeMillis();
		DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Calendar cal = Calendar.getInstance();
		System.out.println("startTime : " + date.format(cal.getTime()));

		Space space = Space3D.getInstance();
		long steps = params.steps;
		boolean isGraphic = params.isGraphic;
        double temperature = params.temperature;
        int numberMolecules = params.numberMolecules;
		double density = params.density;
		double HSDiameter=params.HSDiameter;
		double dipoleStrength = Math.sqrt(params.dipoleStrength2);
		double dielectricOutside = params.dielectricOutside;

		double V = numberMolecules/density;
   		double boxSize = Math.pow(V, 1.0/3.0);
		final double truncation = boxSize * 0.49;
//		System.out.println("******************* dipolar HS, dielectric constant, NVT********************");
		System.out.println("steps = "+steps);
		System.out.println("number of molecules= "+numberMolecules);
		System.out.println("temperature(sim) = "+temperature);
		System.out.println("density(sim)1/angstrom^3 = "+density);
		System.out.println("box size(angstrom): "+boxSize);
		System.out.println("truncation= "+truncation);
//		System.out.println("HSDiameter= "+HSDiameter);
		System.out.println("dipoleStrength = "+dipoleStrength);
		System.out.println("dielectricOutside = "+dielectricOutside);

		final DHSFasterer sim = new DHSFasterer(space,numberMolecules,HSDiameter,dipoleStrength,dielectricOutside,boxSize,temperature,truncation);

		if (isGraphic){
			SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);
			simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
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


		//TODO
		sim.integrator.getMoveManager().setEquilibrating(false);        //set the stepsize
//		System.out.println("equilibration finished");
		//TODO

		int blockNumber = 100;
		int sampleAtInterval = numberMolecules;
        long samplePerBlock = steps / sampleAtInterval / blockNumber;
//		System.out.println("number of blocks is : " + blockNumber);
//		System.out.println("sample per block is : " + samplePerBlock);

		// dipoleSumSquared
        MeterDipoleSumSquared1site dipoleSumSquaredMeter = new MeterDipoleSumSquared1site(space, sim.box, dipoleStrength);
        AccumulatorAverage dipoleSumSquaredAccumulator = new AccumulatorAverageFixed(samplePerBlock);
		DataPumpListener dipolePump = new DataPumpListener(dipoleSumSquaredMeter,dipoleSumSquaredAccumulator, sampleAtInterval);
		sim.integrator.getEventManager().addListener(dipolePump);

		// energy
		MeterPotentialEnergyFromIntegratorFasterer energyMeter = new MeterPotentialEnergyFromIntegratorFasterer(sim.integrator);
		AccumulatorAverage energyAccumulator = new AccumulatorAverageFixed(10);
		DataPumpListener energyPump = new DataPumpListener(energyMeter, energyAccumulator);
		energyAccumulator.setBlockSize(50);
		sim.integrator.getEventManager().addListener(energyPump);

		//AEE
		DipoleSourceDHS dipoleDHS = new DipoleSourceDHS(space, dipoleStrength);// add reaction field potential
		MeterDipoleSumSquaredMappedAverageFasterer AEEMeter = new MeterDipoleSumSquaredMappedAverageFasterer(sim.box, sim.getSpeciesManager(), dipoleStrength, temperature, sim.potentialMaster);
		AEEMeter.setDipoleSource(dipoleDHS);
		AccumulatorAverageCovariance AEEAccumulator = new AccumulatorAverageCovariance(samplePerBlock, true);
		DataPumpListener AEEPump = new DataPumpListener(AEEMeter, AEEAccumulator, sampleAtInterval);

        sim.integrator.getEventManager().addListener(AEEPump);
		sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator, steps));

		//calculate dipoleSumSquared average
        double dipoleSumSquared = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.AVERAGE.index)).x;
        double dipoleSumSquaredERR = ((DataDouble) ((DataGroup) dipoleSumSquaredAccumulator.getData()).getData(dipoleSumSquaredAccumulator.ERROR.index)).x;
        //TODO
        double sum0 = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.AVERAGE.index).getValue(0);
        double ERsum0 = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(0);
        double sum1 = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.AVERAGE.index).getValue(1);
        double ERsum1 = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.ERROR.index).getValue(1);

        IData covariance = ((DataGroup) AEEAccumulator.getData()).getData(AEEAccumulator.BLOCK_COVARIANCE.index);
        covariance.getValue(1);
//		double AEE = sum0 + sum1*sum1;
//		double AEEER = Math.sqrt(ERsum0*ERsum0 + 4*sum1*sum1*ERsum1*ERsum1 -
//				2*ERsum0*sum1*2*ERsum1*covariance.getValue(1)/Math.sqrt(covariance.getValue(0)*covariance.getValue(3)));

        //TODO
		double AEE = sum0;
		double AEEER = ERsum0;
		//TODO


        double volume = sim.box.getBoundary().volume();
		double dipoleFac = 4 * Math.PI * dipoleSumSquared/9.0/volume/temperature;
		double dielectricOutsideFac = dielectricOutside<Double.POSITIVE_INFINITY ? 2*(dielectricOutside-1)/(2*dielectricOutside+1) : 1;
		double x1 =  dipoleSumSquared;
		double B  = dielectricOutsideFac;
		double D  = 4 * Math.PI/9.0/volume/temperature;
//		System.out.println("B:  "+B);
//		System.out.println("D:  "+D);
//		System.out.println("x1:  "+x1);

		double dEpsilondx_1=(B*D+2*D)/(B*D*x1-D*x1+1);
		double dEpsilondx_2=(B*D-D)*(B*D*x1+2D*x1+1)/Math.pow((B*D*x1-D*x1+1),2);
		double dEpsilondx =dEpsilondx_1-dEpsilondx_2;
		double epsilon = (1+B*D*x1+2*D*x1)/(1+B*D*x1-D*x1);
		double epsilonERR_alt = dEpsilondx*dipoleSumSquaredERR;
		double epsilonERR = 3 * D / Math.pow((B*D*x1-D*x1+1),2)*dipoleSumSquaredERR;

//		System.out.println("epsilonERR_alt: "+epsilonERR_alt);

		System.out.println("-<M^2>*bt*bt:\t"+ (-dipoleSumSquared/temperature/temperature) + " -<M^2>*bt*bt_err:\t"+(dipoleSumSquaredERR/temperature/temperature) );
		System.out.println("AEE_new:\t"+ AEE + " AEE_err:\t" + AEEER );
//		System.out.println("Epsilon="+epsilon + "   with err: "+epsilonERR);
//		System.out.println("========================");
		double C = dipoleFac;
		double A = C/(1+B*C);
		double dielectricConstant = (1+2*A)/(1-A);
//		System.out.println("(epsilon-1)/(epsilon+2): "+A);
//		System.out.println("dielectric constant is:  "+dielectricConstant);

		double avgPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.AVERAGE.index)).x / numberMolecules;
		double errPE = ((DataDouble) ((DataGroup) energyAccumulator.getData()).getData(energyAccumulator.ERROR.index)).x / numberMolecules;
//		System.out.println("PE/epsilon: "+avgPE+" "+errPE);
        long endTime = System.currentTimeMillis();
        System.out.println("endTime : " + date.format(cal.getTime()));
        System.out.println("Time taken (in mins): " + (endTime - startTime) / (1000.0 * 60.0));
    }

    public static class DipoleSourceDHS implements DipoleSourceAtomic {//for potential reaction field
        protected final Vector dipoleVector;
        protected double dipoleStrength;

        public DipoleSourceDHS(Space space, double dipole) {
            dipoleStrength = dipole;
            dipoleVector = space.makeVector();
        }

        public Vector getDipole(IAtom atom) {
            IAtomOriented a = (IAtomOriented) atom;
            dipoleVector.E(a.getOrientation().getDirection());
            dipoleVector.TE(dipoleStrength);
            return dipoleVector;
        }
    }

	// ******************* parameters **********************// 
	public static class Param extends ParameterBase {
		public boolean isGraphic = false;
		public double temperature = 1;				
		public int numberMolecules = 100;
		public double density = 0.001;
		public double dipoleStrength2 = 0.25;
		public double HSDiameter = 1.0;
		public double dielectricOutside = 1.0E11;
		public long steps = 1000000;
	}
}
