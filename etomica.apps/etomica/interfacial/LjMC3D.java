/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.interfacial;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.AtomSourceRandomSpecies;
import etomica.atom.DiameterHashByType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageFixed;
import etomica.data.AccumulatorHistory;
import etomica.data.DataDump;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataPumpListener;
import etomica.data.DataSourceCountSteps;
import etomica.data.meter.MeterNMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterPotentialEnergyFromIntegrator;
import etomica.data.meter.MeterProfileByVolume;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.listener.IntegratorListenerAction;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.Potential1;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.HistoryCollapsingAverage;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

/**
 * Simple Lennard-Jones molecular dynamics simulation in 3D
 */
 
public class LjMC3D extends Simulation {
    
    public final PotentialMasterCell potentialMasterCell;
    public final ActivityIntegrate ai;
    public IntegratorMC integrator;
    public SpeciesSpheresMono speciesFluid, speciesTopWall, speciesBottomWall;
    public IBox box;
    public P2SoftSphericalTruncatedForceShifted pFF, pTW, pBW;
    public ConfigurationLammps config;


    public LjMC3D(double temperature, double Lxy, String lammpsFile) {
        super(Space3D.getInstance());
        BoundaryRectangularSlit boundary = new BoundaryRectangularSlit(2, space);
        IVectorMutable bs = (IVectorMutable)boundary.getBoxSize();
        bs.setX(0, Lxy);
        bs.setX(1, Lxy);
        boundary.setBoxSize(bs);
        box = new Box(boundary, space);
        addBox(box);
        
        speciesFluid = new SpeciesSpheresMono(space, new ElementSimple("F"));
        addSpecies(speciesFluid);
        speciesTopWall = new SpeciesSpheresMono(space, new ElementSimple("TW"));
        addSpecies(speciesTopWall);
        speciesBottomWall = new SpeciesSpheresMono(space, new ElementSimple("BW"));
        addSpecies(speciesBottomWall);
        
        config = new ConfigurationLammps(space, lammpsFile, speciesTopWall, speciesBottomWall, speciesFluid);
        config.initializeCoordinates(box);
        
        potentialMasterCell = new PotentialMasterCell(this, 5.49925, space);
        potentialMasterCell.setCellRange(2);
        integrator = new IntegratorMC(this, potentialMasterCell);
        integrator.setTemperature(temperature);
        MCMoveAtom mcMoveAtom = new MCMoveAtom(random, potentialMasterCell, space);
        mcMoveAtom.setAtomSource(new AtomSourceRandomSpecies(getRandom(), speciesFluid));
        MCMoveAtom mcMoveAtomBig = new MCMoveAtom(random, potentialMasterCell, space);
        mcMoveAtomBig.setAtomSource(new AtomSourceRandomSpecies(getRandom(), speciesFluid));
        mcMoveAtomBig.setStepSize(0.5*Lxy);
        ((MCMoveStepTracker)mcMoveAtomBig.getTracker()).setTunable(false);
        
        integrator.getMoveManager().addMCMove(mcMoveAtom);
        integrator.getMoveManager().addMCMove(mcMoveAtomBig);
        
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);

        pFF = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.0, 1.0), 2.5);
        IAtomType leafType = speciesFluid.getLeafType();
        potentialMasterCell.addPotential(pFF,new IAtomType[]{leafType,leafType});

        pBW = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.09985, 0.4), 5.49925);
        potentialMasterCell.addPotential(pBW,new IAtomType[]{leafType,speciesBottomWall.getLeafType()});
        
        pTW = new P2SoftSphericalTruncatedForceShifted(space, new P2LennardJones(space, 1.5, 0.1), 1.68);
        potentialMasterCell.addPotential(pTW,new IAtomType[]{leafType,speciesTopWall.getLeafType()});

        Potential1 p1F = new Potential1(space) {
            
            public double energy(IAtomList atoms) {
                IVector p = atoms.getAtom(0).getPosition();
                double z = p.getX(2);
                double Lz = boundary.getBoxSize().getX(2);
                return (Math.abs(z) > 0.5*Lz) ? Double.POSITIVE_INFINITY : 0;
            }
        };
        potentialMasterCell.addPotential(p1F,new IAtomType[]{leafType});
        
        integrator.setBox(box);

        integrator.getMoveEventManager().addListener(potentialMasterCell.getNbrCellManager(box).makeMCMoveListener());
		
        potentialMasterCell.getNbrCellManager(box).assignCellAll();
    }
    
    public static void main(String[] args) {

        // according to Mastny & de Pablo
        // http://link.aip.org/link/doi/10.1063/1.2753149
        // triple point
        // T = 0.694
        // liquid density = 0.845435
        
        // Agrawal and Kofke:
        //      small      large
        // T    0.698    0.687(4)
        // p    0.0013   0.0011
        // rho  0.854    0.850
        
        // Orkoulas, http://link.aip.org/link/doi/10.1063/1.4758698 says
        // T = 0.7085(5)
        // P = 0.002264(17)
        // rhoL = 0.8405(3)
        // rhoFCC = 0.9587(2)
        // rhoV = 0.002298(18)

        LjMC3DParams params = new LjMC3DParams();
        ParseArgs.doParseArgs(params, args);
        if (args.length==0) {
            params.graphics = false;
            params.lammpsFile = "eq.data.2";
            params.steps = 10000;
            params.T = 0.8;
            params.Lxy = 15.239984;
        }

        final double temperature = params.T;
        final double Lxy = params.Lxy;
        long steps = params.steps;
        boolean graphics = params.graphics;
        String lammpsFile = params.lammpsFile;

        if (!graphics) {
            System.out.println("Running MC with T="+temperature);
            System.out.println(steps+" steps");
        }

        final LjMC3D sim = new LjMC3D(temperature, Lxy, lammpsFile);

        MeterPotentialEnergyFromIntegrator meterPE = new MeterPotentialEnergyFromIntegrator();
        meterPE.setIntegrator(sim.integrator);
        DataFork forkPE = new DataFork();
        DataPumpListener pumpPE = new DataPumpListener(meterPE, forkPE, 10);
        sim.integrator.getEventManager().addListener(pumpPE);

        MeterWallForce meterWF = new MeterWallForce(sim.space, sim.potentialMasterCell, sim.box, sim.speciesTopWall);
        DataFork forkWF = new DataFork();
        DataPumpListener pumpWF = new DataPumpListener(meterWF, forkWF, 1000);
        sim.integrator.getEventManager().addListener(pumpWF);
        
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(sim.potentialMasterCell);
        meterPE2.setBox(sim.box);
        double u = meterPE2.getDataAsScalar();
        System.out.println("Potential energy: "+u);
        System.out.println("Wall force: "+meterWF.getDataAsScalar());
        
        MeterProfileByVolume densityProfileMeter = new MeterProfileByVolume(sim.space);
        densityProfileMeter.setProfileDim(2);
        densityProfileMeter.setBox(sim.box);
        MeterNMolecules meterNMolecules = new MeterNMolecules();
        meterNMolecules.setSpecies(sim.speciesFluid);
        densityProfileMeter.setDataSource(meterNMolecules);
        AccumulatorAverageFixed densityProfileAvg = new AccumulatorAverageFixed(10);
        densityProfileAvg.setPushInterval(10);
        DataPump profilePump = new DataPumpListener(densityProfileMeter, densityProfileAvg, 1000);
        DataDump profileDump = new DataDump();
        densityProfileAvg.addDataSink(profileDump, new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
        IntegratorListenerAction profilePumpListener = new IntegratorListenerAction(profilePump);
        sim.integrator.getEventManager().addListener(profilePumpListener);
        profilePumpListener.setInterval(10);
        
        
        if (graphics) {
            final String APP_NAME = "LjMC3D";
            final SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 3, sim.getSpace(), sim.getController());

            List<DataPump> dataStreamPumps = simGraphic.getController().getDataStreamPumps();
            dataStreamPumps.add(profilePump);

            simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

            simGraphic.makeAndDisplayFrame(APP_NAME);
            DiameterHashByType dh = (DiameterHashByType)simGraphic.getDisplayBox(sim.box).getDiameterHash();
            dh.setDiameter(sim.speciesFluid.getLeafType(), 1.0);
            dh.setDiameter(sim.speciesBottomWall.getLeafType(), 1.09885);
            dh.setDiameter(sim.speciesTopWall.getLeafType(), 1.5);
            
            DataSourceCountSteps dsSteps = new DataSourceCountSteps(sim.integrator);

            AccumulatorHistory historyPE = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyPE.setTimeDataSource(dsSteps);
            forkPE.addDataSink(historyPE);
            DisplayPlot plotPE = new DisplayPlot();
            historyPE.setDataSink(plotPE.getDataSet().makeDataSink());
            plotPE.setLabel("PE");
            simGraphic.add(plotPE);
            
            AccumulatorHistory historyWF = new AccumulatorHistory(new HistoryCollapsingAverage());
            historyWF.setTimeDataSource(dsSteps);
            forkWF.addDataSink(historyWF);
            DisplayPlot plotWF = new DisplayPlot();
            historyWF.setDataSink(plotWF.getDataSet().makeDataSink());
            plotWF.setLabel("Force");
            simGraphic.add(plotWF);
            

            DisplayPlot profilePlot = new DisplayPlot();
            densityProfileAvg.addDataSink(profilePlot.getDataSet().makeDataSink(), new AccumulatorAverage.StatType[]{densityProfileAvg.AVERAGE});
            profilePlot.setLabel("density");
            simGraphic.add(profilePlot);
            
            return;
        }

        long bs = steps/100000;
        if (bs==0) bs=1;
        AccumulatorAverageFixed accPE = new AccumulatorAverageFixed(bs);
        forkPE.addDataSink(accPE);
        AccumulatorAverageFixed accWF = new AccumulatorAverageFixed(bs);
        forkWF.addDataSink(accWF);

        sim.ai.setMaxSteps(steps);
        sim.getController().actionPerformed();
        
        u = meterPE2.getDataAsScalar();
        System.out.println("Potential energy: "+u);
        System.out.println("Wall force: "+meterWF.getDataAsScalar());
        
        double avgPE = accPE.getData().getValue(accPE.AVERAGE.index);
        double errPE = accPE.getData().getValue(accPE.ERROR.index);
        double corPE = accPE.getData().getValue(accPE.BLOCK_CORRELATION.index);
        double avgWF = accWF.getData().getValue(accPE.AVERAGE.index);
        double errWF = accWF.getData().getValue(accPE.ERROR.index);
        double corWF = accWF.getData().getValue(accPE.BLOCK_CORRELATION.index);
        
        if (steps>100000) {
            System.out.println(String.format("Average potential energy: %25.15e %10.4e % 5.3f\n",avgPE,errPE,corPE));
            System.out.println(String.format("Average wall force: %25.15e %10.4e % 5.3f\n",avgWF,errWF,corWF));
        }
        else {
            System.out.println("Average potential energy: "+avgPE);
            System.out.println("Average wall force: "+avgWF);
        }
        
        WriteConfigurationInterfacial configWriter = new WriteConfigurationInterfacial(sim.space);
        configWriter.setSpecies(sim.speciesFluid);
        IVectorMutable unshift = sim.space.makeVector();
        unshift.Ea1Tv1(-1, sim.config.getShift());
        configWriter.setShift(unshift);
        configWriter.setBox(sim.box);
        configWriter.setFileName("xyz_000.dat");
        configWriter.actionPerformed();
    }
    
    public static class ConfigurationLammps implements Configuration {
        protected final ISpace space;
        protected final String filename;
        protected final ISpecies[] species;
        protected IVectorMutable shift;
        
        public ConfigurationLammps(ISpace space, String filename, ISpecies topWall, ISpecies bottomWall, ISpecies fluid) {
            this.space = space;
            this.filename = filename;
            this.species = new ISpecies[]{topWall, fluid, bottomWall};
            shift = space.makeVector();
        }
        
        public void initializeCoordinates(IBox box) {
            List<IVectorMutable>[] coords = new List[3];
            coords[0] = new ArrayList<IVectorMutable>();
            coords[1] = new ArrayList<IVectorMutable>();
            coords[2] = new ArrayList<IVectorMutable>();
            double zMin = Double.POSITIVE_INFINITY;
            double zMax = -Double.POSITIVE_INFINITY;
            try {
                FileReader fr = new FileReader(filename);
                BufferedReader bufReader = new BufferedReader(fr);
                String line = null;
                while ((line = bufReader.readLine()) != null) {
                    String[] bits = line.split("\\t");
                    if (bits.length != 5) continue;
                    int aType = Integer.parseInt(bits[1]);
                    IVectorMutable xyz = space.makeVector();
                    for (int i=0; i<3; i++) {
                        xyz.setX(i, Double.parseDouble(bits[2+i]));
                    }
                    if (xyz.getX(2) > zMax) zMax = xyz.getX(2);
                    if (xyz.getX(2) < zMin) zMin = xyz.getX(2);
                    coords[aType-1].add(xyz);
                }
                bufReader.close();
            }
            catch (IOException ex) {
                throw new RuntimeException(ex);
            }
            double Lz = zMax-zMin + 0.8;
            IVectorMutable bs = (IVectorMutable)box.getBoundary().getBoxSize();
            bs.setX(2, Lz);
            box.getBoundary().setBoxSize(bs);
            shift.setX(0, -0.5*bs.getX(0));
            shift.setX(1, -0.5*bs.getX(1));
            shift.setX(2, -0.5*Lz - zMin);
            for (int i=0; i<3; i++) {
                box.setNMolecules(species[i], coords[i].size());
                IMoleculeList m = box.getMoleculeList(species[i]);
                for (int j=0; j<coords[i].size(); j++) {
                    IVectorMutable p = m.getMolecule(j).getChildList().getAtom(0).getPosition();
                    p.Ev1Pv2(coords[i].get(j), shift);
                }
            }
            
        }
        
        public IVector getShift() {
            return shift;
        }
    }

    public static class LjMC3DParams extends ParameterBase {
        public double T = 2.0;
        public long steps = 100000;
        public boolean graphics = false;
        public double Lxy = 0;
        public String lammpsFile = "";
    }
}
