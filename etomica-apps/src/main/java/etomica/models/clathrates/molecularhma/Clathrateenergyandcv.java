
package etomica.models.clathrates.molecularhma;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationFileBinary;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCovariance;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.DataPumpListener;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergysquare;
import etomica.exception.ConfigurationOverlapException;
import etomica.integrator.IntegratorVelocityVerletRattle;
import etomica.integrator.IntegratorVelocityVerletShake;
import etomica.models.clathrates.molecularhma.MinimizationTIP4P.ChargeAgentSourceRPM;
import etomica.models.water.ConformationWaterTIP4P;
import etomica.models.water.SpeciesWater4P;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculePositionCOM;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.units.Calorie;
import etomica.units.Mole;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Clathrateenergyandcv extends Simulation {
    protected static double[] initialU;
    protected Box box;
    protected PotentialMaster potentialMaster;
    protected SpeciesWater4P species;
    protected MeterPotentialEnergy meterPE;
    protected Potential2SoftSpherical potentialLJ;
    protected Potential2SoftSphericalLS potentialLJLS;
    protected EwaldSummation potentialES;
    public IntegratorVelocityVerletRattle integrator;
    //   protected final AtomLeafAgentManager<Vector> forceManager;
    public ActivityIntegrate ai;
    protected final MoleculeAgentManager latticeCoordinates;

    public Clathrateenergyandcv(Space space, double temperature, int numCells, double rCutRealES, double rCutLJ, boolean isIce, double kCut, double shakeTol, boolean unitCells, final boolean doTranslation, final boolean doRotation, double timeInterval) {

        super(space);
        box = new Box(space);
        species = new SpeciesWater4P(getSpace(), true);
        addSpecies(species);
        addBox(box);

        box.setNMolecules(species, 46 * numCells * numCells * numCells);
        box.setDensity(46 / 12.03 / 12.03 / 12.03);
        ChargeAgentSourceRPM agentSource = new ChargeAgentSourceRPM(species, isIce);
        AtomLeafAgentManager<EwaldSummation.MyCharge> atomAgentManager = new AtomLeafAgentManager<EwaldSummation.MyCharge>(agentSource, box);
        int[] nC = new int[]{numCells, numCells, numCells};
        if (unitCells) {
            numCells = 1;
        }
        ConfigurationFile config = new ConfigurationFile(numCells + "ncfinalPos");
//////////minimized file
        if (unitCells) {
            ConfigurationFileBinary.replicate(config, box, nC, space);
        } else {
            config.initializeCoordinates(box);
        }
        double a0 = box.getBoundary().getBoxSize().getX(0);
        double[] rC = new double[]{a0, a0, a0};
        double sigma, epsilon; //TIP4P
        if (isIce) {
            sigma = 3.1668;
            epsilon = 106.1;//TIP4P/Ice
        } else {//TIP4P
            double A = 600E3; // kcal A^12 / mol
            double C = 610.0; // kcal A^6 / mol
            double s6 = A / C;
            sigma = Math.pow(s6, 1.0 / 6.0);
            epsilon = Mole.UNIT.fromSim(Calorie.UNIT.toSim(C / s6 * 1000)) / 4.0;
        }

        latticeCoordinates = new MoleculeAgentManager(this, box, new MoleculeSiteSource(space, new MoleculePositionCOM(space), new Clathrateenergyandcv.WaterOrientationDefinition(space)));
 //       AtomLeafAgentManager<Vector> atomLatticeCoordinates = new AtomLeafAgentManager<>(new AtomSiteSource(space), box, Vector.class);

         potentialES = new EwaldSummation(box, atomAgentManager, space, kCut, rCutRealES);

        potentialLJ = new P2LennardJones(space, sigma, epsilon);
        potentialLJLS = new Potential2SoftSphericalLS(space, rCutLJ, rC, potentialLJ);
         potentialMaster = new PotentialMaster();
        potentialMaster.addPotential(potentialES, new AtomType[0]);
        potentialMaster.addPotential(potentialLJLS, new AtomType[]{species.getOxygenType(), species.getOxygenType()});
        int maxIterations = 100;
        integrator = new IntegratorVelocityVerletRattle(this, potentialMaster, box);
        integrator.setShakeTolerance(shakeTol);
        double lOH = ConformationWaterTIP4P.bondLengthOH;
        double lHH = Math.sqrt(2 * lOH * lOH * (1 - Math.cos(ConformationWaterTIP4P.angleHOH)));
        double lOM = ConformationWaterTIP4P.rOM;
        double lMH = Math.sqrt(lOH * lOH + lOM * lOM - 2 * lOH * lOM * Math.cos(0.5 * ConformationWaterTIP4P.angleHOH));
        IntegratorVelocityVerletShake.BondConstraints bondConstraints = new IntegratorVelocityVerletShake.BondConstraints(new int[][]{{0, 2}, {1, 2}, {0, 1}}, new double[]{lOH, lOH, lHH}) {

            Vector vectorSum = space.makeVector();
            Vector oVector = space.makeVector();
            Vector centerMass = space.makeVector();
            Vector h1Vector = space.makeVector();
            Vector h2Vector = space.makeVector();
            Vector mVector = space.makeVector();
            //            Vector newForce = space.makeVector();
//            Vector newTorque = space.makeVector();
            Vector dr = space.makeVector();
//make force on m site=0 as its mass is infinite so acceleratn is infinite, but keep force,torque of com same
            public void redistributeForces(IMolecule molecule, AtomLeafAgentManager<Vector> agentManager) {
                IAtomList leafList = molecule.getChildList();
                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                Vector h1torque = space.makeVector();
                Vector h2torque = space.makeVector();
                Vector h1h2 = space.makeVector();
                Vector om = space.makeVector();
                Vector oTorque = space.makeVector();
                Vector mTorque = space.makeVector();
                Vector totalTorque = space.makeVector();
                Vector totalForce = space.makeVector();
                double hMass = leafList.get(0).getType().getMass();
                double oMass = leafList.get(2).getType().getMass();
                double hMassPercent = hMass / (2 * hMass + oMass);
                double oMassPercent = oMass / (2 * hMass + oMass);
                centerMass.Ea1Tv1(hMass, h1);
                centerMass.PEa1Tv1(hMass, h2);
                centerMass.PEa1Tv1(oMass, o);
                centerMass.TE(1 / (2 * hMass + oMass));
                mVector.Ev1Mv2(m, centerMass);
                h1Vector.Ev1Mv2(h1, centerMass);
                h2Vector.Ev1Mv2(h2, centerMass);
                oVector.Ev1Mv2(o, centerMass);
                vectorSum.E(h1Vector);
                vectorSum.PE(h2Vector);
                vectorSum.PEa1Tv1(-2.0, oVector);
                Vector h1Force = agentManager.getAgent(leafList.get(0));
                Vector h2Force = agentManager.getAgent(leafList.get(1));
                Vector oForce = agentManager.getAgent(leafList.get(2));
                Vector mForce = agentManager.getAgent(leafList.get(3));

                boolean testForceAndTorque = false;
                Vector totalTorqueNew = space.makeVector();
                Vector totalForceNew = space.makeVector();
                if (testForceAndTorque) {
                    totalForce.E(h1Force);
                    totalForce.PE(h2Force);
                    totalForce.PE(mForce);
                    totalForce.PE(oForce);// test for total force

                    h1torque.E(h1Force);
                    h1torque.XE(h1Vector);
                    h2torque.E(h2Force);
                    h2torque.XE(h2Vector);
                    oTorque.E(oForce);
                    oTorque.XE(oVector);
                    mTorque.E(mForce);
                    mTorque.XE(mVector);
                    totalTorque.E(h1torque);
                    totalTorque.PE(h2torque);
                    totalTorque.PE(oTorque);
                    totalTorque.PE(mTorque);//test for total torque
                }

                h1h2.Ev1Mv2(h2, h1);
                om.Ev1Mv2(m, o);
                Vector a0 = space.makeVector();
                Vector a1 = space.makeVector();
                Vector a2 = space.makeVector();
                a0.E(mVector);
                a0.normalize();
                a1.Ea1Tv1(-1, h1h2);   //#########WHAT DOES Ea1Tv1 DO?????
                a1.normalize();
                a2.E(a0);//equal
                a2.XE(a1);//cross a2=a2Xa1

                //separate mForce into 3 direction
                Vector f0 = space.makeVector();
                Vector f1 = space.makeVector();
                Vector f2 = space.makeVector();
                f0.Ea1Tv1(mForce.dot(a0), a0);
                f1.Ea1Tv1(mForce.dot(a1), a1);
                f2.Ea1Tv1(mForce.dot(a2), a2);

                //Translation part
                h1Force.PEa1Tv1(hMassPercent, mForce);
                h2Force.PEa1Tv1(hMassPercent, mForce);
                oForce.PEa1Tv1(oMassPercent, mForce);

                //Rotation part
                dr.Ea1Tv1(f1.dot(a1) * mVector.dot(a0) / (-2 * oVector.dot(a0) + 2 * h1Vector.dot(a0)), a1);
//                System.out.println("part1="+ f1.dot(mVector)/(-2*oVector.dot(a0)+2*h1Vector.dot(a0)));
//                System.out.println("part1Force={"+dr.getX(0)+","+dr.getX(1)+","+dr.getX(2)+"};" );
                h1Force.PE(dr);
                h2Force.PE(dr);
                oForce.PEa1Tv1(-2, dr);

                dr.Ea1Tv1(f2.dot(a2) * mVector.dot(a0) / (-2 * oVector.dot(a0) + 2 * h1Vector.dot(a0)), a2);
//                System.out.println("part2="+f2.dot(mVector)/(-2*oVector.dot(a0)+2*h1Vector.dot(a0)));
//                System.out.println("part2Force={"+dr.getX(0)+","+dr.getX(1)+","+dr.getX(2)+"};" );
                h1Force.PE(dr);
                h2Force.PE(dr);
                oForce.PEa1Tv1(-2, dr);

                mForce.E(0);//remove mForce at the end

                if (testForceAndTorque) {
                    //test for total force
                    totalForceNew.E(h1Force);
                    totalForceNew.PE(h2Force);
                    totalForceNew.PE(mForce);
                    totalForceNew.PE(oForce);//New total force

                    System.out.println("totalForceOld ={ " + totalForce.getX(0) + "," + totalForce.getX(1) + "," + totalForce.getX(2) + "};");
                    System.out.println("totalForceNew ={ " + totalForceNew.getX(0) + "," + totalForceNew.getX(1) + "," + totalForceNew.getX(2) + "};");


                    //test for total torque
                    h1torque.E(h1Force);
                    h1torque.XE(h1Vector);
                    h2torque.E(h2Force);
                    h2torque.XE(h2Vector);
                    oTorque.E(oForce);
                    oTorque.XE(oVector);
                    mTorque.E(mForce);
                    mTorque.XE(mVector);
                    totalTorqueNew.E(h1torque);
                    totalTorqueNew.PE(h2torque);
                    totalTorqueNew.PE(mTorque);
                    totalTorqueNew.PE(oTorque);//New total torque
                    System.out.println("totalTorqueOld ={ " + totalTorque.getX(0) + "," + totalTorque.getX(1) + "," + totalTorque.getX(2) + "};");
                    System.out.println("totalTorqueNew ={ " + totalTorqueNew.getX(0) + "," + totalTorqueNew.getX(1) + "," + totalTorqueNew.getX(2) + "};");
                    System.exit(2);
                }


                //redistribute the forces s.t. only translation
                if (!doRotation) {
                    totalForce.E(h1Force);
                    totalForce.PE(h2Force);
                    totalForce.PE(mForce);
                    totalForce.PE(oForce);
                    h1Force.E(0);
                    h2Force.E(0);
                    oForce.E(0);
                    mForce.E(0);
                    h1Force.PEa1Tv1(hMassPercent, totalForce);
                    h2Force.PEa1Tv1(hMassPercent, totalForce);
                    oForce.PEa1Tv1(oMassPercent, totalForce);
                }


                //redistribute s.t. only rotation
                if (!doTranslation) {
                    totalForce.E(h1Force);
                    totalForce.PE(h2Force);
                    totalForce.PE(mForce);
                    totalForce.PE(oForce);
                    h1Force.PEa1Tv1(-1.0 * hMassPercent, totalForce);
                    h2Force.PEa1Tv1(-1.0 * hMassPercent, totalForce);
                    oForce.PEa1Tv1(-1.0 * oMassPercent, totalForce);
                }

            }

            public void relaxMolecule(IMolecule molecule) {
                IAtomList leafList = molecule.getChildList();

                Vector h1 = leafList.get(0).getPosition();
                Vector h2 = leafList.get(1).getPosition();
                Vector o = leafList.get(2).getPosition();
                Vector m = leafList.get(3).getPosition();
                m.E(h1);
                m.PE(h2);
                m.PEa1Tv1(-2, o);
                m.TE(0.15 / Math.sqrt(m.squared())); // TIP4P
//            	m.TE(0.1577/Math.sqrt(m.squared())); TIP4P/Ice
                m.PE(o);

            }

        };
        integrator.getShakeAgentManager().setAgent(species, bondConstraints);
        integrator.setTimeStep(timeInterval);
        integrator.setMaxIterations(maxIterations);
//        integrator.setOrientAtom((IAtom)((IMolecule)box.getMoleculeList(speciesOrient).getAtom(0)).getChildList().getAtom(0));
        integrator.setIsothermal(true);
        integrator.setTemperature(temperature);

        // WHAT'S NEED TO CONVERT UNITS

        integrator.setThermostatInterval(100);
        integrator.setThermostatNoDrift(true);

        try {
            integrator.reset();
        } catch (ConfigurationOverlapException e) {
        }
        ai = new ActivityIntegrate(integrator);
//        System.out.println("using rigid with dt="+dt);
        getController().addAction(ai);
//        System.out.println("h1 at "+((IAtomPositioned)box.getLeafList().getAtom(0)).getPosition());
//        System.out.println("o at "+((IAtomPositioned)box.getLeafList().getAtom(2)).getPosition());

    }


    public static void main(String[] args) {


        final long startTime = System.currentTimeMillis();
        SimParams params = new SimParams();
        ParseArgs.doParseArgs(params, args);
        final double temperature = params.temperature;
        int numCells = params.numCells;
        int numSteps = params.numSteps;
        double rCutRealES = params.rCutRealES;
        double rCutLJ = params.rCutLJ;
        boolean isIce = params.isIce;
        double kCut = params.kCut;
        double shakeTol = params.shakeTol;
        boolean uniteCells = params.unitCells;
        boolean doTranslation = params.doTranslation;
        boolean doRotation = params.doRotation;
         double timeInterval = params.timeInterval;
        boolean doMapping = params.doMapping;
        boolean doConventional = params.doConventional;
        final Clathrateenergyandcv sim = new Clathrateenergyandcv(Space3D.getInstance(), temperature, numCells, rCutRealES, rCutLJ, isIce, kCut, shakeTol, uniteCells, doTranslation, doRotation, timeInterval);
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE2.setBox(sim.box);
        final double latticeEnergy = meterPE2.getDataAsScalar();
        System.out.println("latticeEnergy = " + latticeEnergy);

        final int[] nC = params.nC;
        final double[] a0 = params.a0;
        int blockSize = numSteps >= 1000 ? (numSteps / 1000) : 1;

        final MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        final MeterPotentialEnergysquare meterPEsquare = new MeterPotentialEnergysquare(sim.potentialMaster);
        meterPEsquare.setBox(sim.box);

        AccumulatorAverageCovariance accumulatorAverageFixedDADB = null;
        DataPumpListener dataPumpListenerDADB = null;
        if (doMapping) {
            MeterDADBWaterTIP4P  meterDADB = new MeterDADBWaterTIP4P(nC,a0, rCutLJ,sim.potentialLJ, sim.potentialLJLS,sim.potentialES,sim.space,meterPE, sim.box, params.nBasis,sim.potentialMaster, temperature, sim.latticeCoordinates);

                  meterDADB.doTranslation = doTranslation;
            meterDADB.doRotation = doRotation;
   //         meterDADB.getData();
   //         System.exit(0);
            accumulatorAverageFixedDADB = new AccumulatorAverageCovariance(blockSize);
            dataPumpListenerDADB = new DataPumpListener(meterDADB, accumulatorAverageFixedDADB, 10);
            //        MeterDADB.justU = true;
        }
        //TODO try lower interval 5-10

        AccumulatorAverageFixed accumulatorAverageFixedPE = null;
        AccumulatorAverageFixed accumulatorAverageFixedPEsquare = null;

        DataPumpListener dataPumpListenerPE = null;
        DataPumpListener dataPumpListenerPEsquare = null;

        if (doConventional) {
            accumulatorAverageFixedPE = new AccumulatorAverageFixed(blockSize);
            dataPumpListenerPE = new DataPumpListener(meterPE, accumulatorAverageFixedPE, 10);
            accumulatorAverageFixedPEsquare = new AccumulatorAverageFixed(blockSize);
            dataPumpListenerPEsquare = new DataPumpListener(meterPEsquare, accumulatorAverageFixedPEsquare, 10);
        }

 //      sim.ai.setMaxSteps(numSteps / 10);
 //       sim.getController().actionPerformed();

        sim.ai.setMaxSteps(numSteps);
        if (doMapping) sim.integrator.getEventManager().addListener(dataPumpListenerDADB);
        if (doConventional) sim.integrator.getEventManager().addListener(dataPumpListenerPE);
        if (doConventional) sim.integrator.getEventManager().addListener(dataPumpListenerPEsquare);

        sim.getController().reset();
        sim.getController().actionPerformed();

        System.out.println("numSteps= " + numSteps);
        System.out.println("temperature= " + temperature);
        IMoleculeList molecules = sim.box.getMoleculeList();
        System.out.println("beginLE = " + latticeEnergy);
        int N = sim.box.getMoleculeList().size();
        double fac = (doRotation ? 1.5 : 0) * N + (doTranslation ? 1.5 : 0) * (N - 1);
        System.out.println("harmonicE= " + fac * (molecules.size() * temperature));
        System.out.println("main:doTranslation:   " + doTranslation + "  doRotation:  " + doRotation);
        System.out.println("timeInterval= " + timeInterval);


        long endTime = System.currentTimeMillis();
        DateFormat date = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
        Calendar cal = Calendar.getInstance();
        System.out.println(date.format(cal.getTime()));
        double totalTime = (endTime - startTime) / (1000.0 * 60.0);
        System.out.println("totalTime= " + totalTime + " mins");

        if (doMapping) {
            double mappingAverageanh = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(0);
            double mappingErroranh = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.ERROR).getValue(0);
            double mappingCoranh = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_CORRELATION).getValue(0);
            double mappingStdevanh = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.STANDARD_DEVIATION).getValue(0);

            double mappingAverageanhsquare = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(1);
            double mappingErroranhsquare = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.ERROR).getValue(1);
            double mappingCoranhsquare = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_CORRELATION).getValue(1);
            double mappingStdevanhsquare = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.STANDARD_DEVIATION).getValue(1);

            double mappingAveragedelrphidelr = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(2);
            double mappingErrordelrphidelr = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.ERROR).getValue(2);
            double mappingCordelrphidelr = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_CORRELATION).getValue(2);
            double mappingStdevdelrphidelr = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.STANDARD_DEVIATION).getValue(2);

            double conPEAverage = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(3);
            double conPEError = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.ERROR).getValue(3);
            double conPECor = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_CORRELATION).getValue(3);
            double conPEStdev = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.STANDARD_DEVIATION).getValue(3);

            double conPEsquareAverage = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(4);
            double conPEsquareError = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.ERROR).getValue(4);
            double conPEsquareCor = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_CORRELATION).getValue(4);
            double conPEsquareStdev = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.STANDARD_DEVIATION).getValue(4);


            double forcorr1 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(5);
            double forcorr2 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(6);
            double forcorr3 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(7);
            double forcorr4 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(8);
            double forcorr5 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(9);
            double forcorr6 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(10);
            double forcorr7 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.AVERAGE).getValue(11);



            //           double mappingcovmatrix0 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(0);
 //           double mappingcovmatrix1 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(1);
 //           double mappingcovmatrix2 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(2);
 //           double mappingcovmatrix3 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(3);
 //           double mappingcovmatrix4 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(4);
  //          double mappingcovmatrix5 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(5);
   //         double mappingcovmatrix6 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(6);
   //         double mappingcovmatrix7 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(7);
  //          double mappingcovmatrix8 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(8);
 //           double mappingcovmatrix9 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(9);
 //           double mappingcovmatrix10 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(10);
 //           double mappingcovmatrix11 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(11);
 //           double mappingcovmatrix12 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(12);
 //           double mappingcovmatrix13 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(13);
  //          double mappingcovmatrix14 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(14);
  //          double mappingcovmatrix15 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(15);
  //          double mappingcovmatrix16 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(16);
  //          double mappingcovmatrix17 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(17);
  //          double mappingcovmatrix18 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(18);
  //          double mappingcovmatrix19 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(19);
  //          double mappingcovmatrix20 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(20);

 //           double mappingcovmatrix21 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(21);
 //           double mappingcovmatrix22 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(22);
 //           double mappingcovmatrix23 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(23);
 //           double mappingcovmatrix24 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(24);
 //           double mappingcovmatrix25 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(25);
 //           double mappingcovmatrix26 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(26);
 //           double mappingcovmatrix27 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(27);
 //           double mappingcovmatrix28 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(28);
            //          double mappingcovmatrix29 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(29);
  //          double mappingcovmatrix30 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(30);

  //          double mappingcovmatrix31 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(31);
  //          double mappingcovmatrix32 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(32);
  //          double mappingcovmatrix33 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(33);
  //          double mappingcovmatrix34 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(34);
  //          double mappingcovmatrix35 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(35);
  //          double mappingcovmatrix36 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(36);
  //          double mappingcovmatrix37 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(37);
  //          double mappingcovmatrix38 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(38);
  //          double mappingcovmatrix39 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(39);
  //          double mappingcovmatrix40 = accumulatorAverageFixedDADB.getData(accumulatorAverageFixedDADB.BLOCK_COVARIANCE).getValue(40);

             System.out.println(" mappingAverageanh=\t" + mappingAverageanh + " mappingErroranh=\t" + mappingErroranh + " mappingCoranh=\t" + mappingCoranh + " mappingstdevanh=\t" + mappingStdevanh + " mappingAverageanhsquare=\t" + mappingAverageanhsquare + " mappingErroranhsquare=\t" + mappingErroranhsquare + " mappingCoranhsquare=\t" + mappingCoranhsquare + " mappingstdevanhsquare=\t" + mappingStdevanhsquare);
            System.out.println(" conPEAverage=\t" + conPEAverage + " conPEError=\t" + conPEError + " conPECor=\t" + conPECor + " conPEStdev=\t" + conPEStdev + " conPEsquareAverage=\t" + conPEsquareAverage + " conPEsquareError=\t" + conPEsquareError + " conPEsquareCor=\t" + conPEsquareCor + " conPEsquareStdev=\t" + conPEsquareStdev );
            System.out.println(" mappingAveragedelrphidelr=\t" + mappingAveragedelrphidelr + " mappingerrordelrphidelr=\t" + mappingErrordelrphidelr + " mappingCordelrphidelr=\t" + mappingCordelrphidelr + " mappingStdevdelrphidelr=\t" +mappingStdevdelrphidelr);

            System.out.println(" forcorr1=\t" +forcorr1+" forcorr2=\t" +forcorr2+" forcorr3=\t" +forcorr3+" forcorr4=\t" +forcorr4+" forcorr5=\t" +forcorr5+" forcorr6=\t" +forcorr6+" forcorr7=\t" +forcorr7);
                    //          System.out.println(" mappingcovmatrix0=" +mappingcovmatrix0+" mappingcovmatrix1=" +mappingcovmatrix1+" mappingcovmatrix2=" +mappingcovmatrix2+" mappingcovmatrix3=" +mappingcovmatrix3+" mappingcovmatrix4=" +mappingcovmatrix4+" mappingcovmatrix5=" +mappingcovmatrix5+" mappingcovmatrix6=" +mappingcovmatrix6+" mappingcovmatrix7=" +mappingcovmatrix7+" mappingcovmatrix8=" +mappingcovmatrix8+" mappingcovmatrix9=" +mappingcovmatrix9+" mappingcovmatrix10=" +mappingcovmatrix10);
  //          System.out.println(" mappingcovmatrix11=" +mappingcovmatrix11+" mappingcovmatrix12=" +mappingcovmatrix12+" mappingcovmatrix13=" +mappingcovmatrix13+" mappingcovmatrix14=" +mappingcovmatrix14+" mappingcovmatrix15=" +mappingcovmatrix15+" mappingcovmatrix16=" +mappingcovmatrix16+" mappingcovmatrix17=" +mappingcovmatrix17+" mappingcovmatrix18=" +mappingcovmatrix18+" mappingcovmatrix19=" +mappingcovmatrix19+" mappingcovmatrix20=" +mappingcovmatrix20);
  //          System.out.println(" mappingcovmatrix21=" +mappingcovmatrix21+" mappingcovmatrix22=" +mappingcovmatrix22+" mappingcovmatrix23=" +mappingcovmatrix23+" mappingcovmatrix24=" +mappingcovmatrix24+" mappingcovmatrix25=" +mappingcovmatrix25+" mappingcovmatrix26=" +mappingcovmatrix26+" mappingcovmatrix27=" +mappingcovmatrix27+" mappingcovmatrix28=" +mappingcovmatrix28+" mappingcovmatrix29=" +mappingcovmatrix29+" mappingcovmatrix30=" +mappingcovmatrix30);
  //          System.out.println(" mappingcovmatrix31=" +mappingcovmatrix31+" mappingcovmatrix32=" +mappingcovmatrix32+" mappingcovmatrix33=" +mappingcovmatrix33+" mappingcovmatrix34=" +mappingcovmatrix34+" mappingcovmatrix35=" +mappingcovmatrix35+" mappingcovmatrix36=" +mappingcovmatrix36+" mappingcovmatrix37=" +mappingcovmatrix37+" mappingcovmatrix38=" +mappingcovmatrix38+" mappingcovmatrix39=" +mappingcovmatrix39+" mappingcovmatrix40=" +mappingcovmatrix40);

        }

        if (doConventional) {
            double PEAverage = accumulatorAverageFixedPE.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double PEAError = accumulatorAverageFixedPE.getData(AccumulatorAverage.ERROR).getValue(0);
            double PEstdev = accumulatorAverageFixedPE.getData(AccumulatorAverage.STANDARD_DEVIATION).getValue(0);
            double PECor = accumulatorAverageFixedPE.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);



            double PEAveragesquare = accumulatorAverageFixedPEsquare.getData(AccumulatorAverage.AVERAGE).getValue(0);
            double PEAErrorsquare = accumulatorAverageFixedPEsquare.getData(AccumulatorAverage.ERROR).getValue(0);
            double PEstdevsquare = accumulatorAverageFixedPEsquare.getData(AccumulatorAverage.STANDARD_DEVIATION).getValue(0);
            double PECorsquare = accumulatorAverageFixedPEsquare.getData(AccumulatorAverage.BLOCK_CORRELATION).getValue(0);

            System.out.println("PEAverage=\t" +PEAverage+" PEerror=\t" +PEAError+ " PECor=\t" + PECor + " PEstdev=\t" + PEstdev + " anharmoniccon "+ (PEAverage - latticeEnergy - fac * temperature) + " time= " + totalTime + " PEAveragesquare=\t" +PEAveragesquare+" PEerrorsquare=\t" +PEAErrorsquare+ " PECorsquare=\t" + PECorsquare + " PEstdevsquare=\t" + PEstdevsquare );
        }


     }


    public static class WaterOrientationDefinition implements MoleculeSiteSource.MoleculeOrientationDefinition {
        protected final OrientationFull3D or;
        protected final Vector v1, v2;

        public WaterOrientationDefinition(Space space) {
            or = new OrientationFull3D(space);
            v1 = space.makeVector();
            v2 = space.makeVector();

        }

        public IOrientation getOrientation(IMolecule molecule) {
            IAtomList leafList = molecule.getChildList();
            Vector h1 = leafList.get(0).getPosition();
            Vector h2 = leafList.get(1).getPosition();
            Vector o = leafList.get(2).getPosition();
            Vector m = leafList.get(3).getPosition();
            v1.Ev1Mv2(m, o);
            v1.normalize();
            v2.Ev1Mv2(h2, h1);
            v2.normalize();
            or.setDirections(v1, v2);
            //v1 is a0 and v2 is a1
            return or;
        }
    }

    public static class SimParams extends ParameterBase {
        //    	public String configFile = "config_from_paper_HHO_shiftedL_2_sI";
        public String configFile = "finalPos";
          public double[] a0 = new double[]{12.03, 12.03, 12.03};//sI
        public int nBasis = 46;//sI
          public boolean includeM = true;
        int nX = 1;
        public int[] nC = new int[]{nX, nX, nX};
        public int numCells = 1;
        public int numSteps = 10;
        public double timeInterval = 0.000000075;
        public double temperature = 50;
        public double rCutLJ = 11;
        public double rCutRealES = 11;
        public double kCut = 1.5;
        public boolean isIce = false;
        public double shakeTol = 1e-12;
         public boolean unitCells = false;
        public boolean doRotation = true;
        public boolean doTranslation = true;
        public boolean doMapping = true;
        public boolean doConventional = true;
    }
}
