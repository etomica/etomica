/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.*;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterVolume;
import etomica.exception.ConfigurationOverlapException;
import etomica.graphics.*;
import etomica.integrator.IntegratorListenerAction;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMove;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.integrator.mcmove.MCMoveVolume;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveOrthorhombicHexagonal;
import etomica.modifier.Modifier;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2HardSphere;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space2d.Space2D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;
import etomica.util.ParameterBase;
import etomica.util.ParseArgs;

import java.awt.*;
import java.io.FileWriter;
import java.io.IOException;

/**
 * NPT simulation for hard sphere solid using an MCMove that does coordinate
 * scaling with volume changes.
 */
public class HSNPT extends Simulation {
    
    public final PotentialMasterList potentialMaster;
    public final IntegratorMC integrator;
    public final SpeciesSpheresMono species;
    public final Box box;
    public final CoordinateDefinition coordinateDefinition;
    public final P2HardSphere pCross;

    public HSNPT(Space _space, int numAtoms, double rho, boolean nvt, boolean fancyMove, double sigma2) {
        super(_space);

        species = new SpeciesSpheresMono(this, space);
        species.setIsDynamic(true);
        addSpecies(species);

        double sigma = 1.0;
        SpeciesSpheresMono species2 = null;
        if (sigma2 != sigma) {
            species2 = new SpeciesSpheresMono(this, space);
            species.setIsDynamic(true);
            addSpecies(species2);
        }

        potentialMaster = new PotentialMasterList(this, space);

        double neighborRangeFac = 1.4;
        double l = Math.pow(numAtoms / rho, 1.0 / 3.0);
        if (_space.D() == 2) {
            int nx = 10;
            int ny = 6;
            if (numAtoms != nx * ny * 2) throw new RuntimeException("oops");
            double bx = 1;
            double v1 = Math.sqrt(3) / rho;
            double v2 = 2 / rho;
            bx *= nx * Math.sqrt(v2 / v1);
            l = bx;
        }
        potentialMaster.setCellRange(1);
        potentialMaster.setRange(neighborRangeFac * sigma);
        box = this.makeBox();
        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0, box);
        this.getController().addActivity(new ActivityIntegrate(integrator));
        AtomType type1 = species.getLeafType();

        P2HardSphere p2 = new P2HardSphere(space, sigma, false);
        potentialMaster.addPotential(p2, new AtomType[]{type1, type1});
        box.setNMolecules(species, numAtoms);
        if (sigma2 != sigma) {
            AtomType type2 = species2.getLeafType();
            pCross = new P2HardSphere(space, (sigma + sigma2) / 2.0, false);
            potentialMaster.addPotential(pCross, new AtomType[]{type1, type2});
            box.setNMolecules(species, numAtoms - 1);
            box.setNMolecules(species2, 1);
        } else {
            pCross = null;
        }

        if (_space.D() == 3) {
            box.getBoundary().setBoxSize(Vector.of(new double[]{l, l, l}));
        } else {
            box.getBoundary().setBoxSize(Vector.of(new double[]{l, l * Math.sqrt(3) * 6 / 10}));
        }
        if (_space.D() == 3) {
            int n = (int) Math.round(Math.pow(numAtoms / 4, 1.0 / 3.0));
            coordinateDefinition = new CoordinateDefinitionLeaf(box, new PrimitiveCubic(space, l / n), new BasisCubicFcc(), space);
            coordinateDefinition.initializeCoordinates(new int[]{n, n, n});
        } else {
            int nx = 10;
            int ny = 6;
            coordinateDefinition = new CoordinateDefinitionLeaf(box, new PrimitiveOrthorhombicHexagonal(space, l / nx), new BasisOrthorhombicHexagonal(), space);
            coordinateDefinition.initializeCoordinates(new int[]{nx, ny});
        }

        potentialMaster.getNeighborManager(box).reset();

        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(potentialMaster, box);
        MCMoveAtomCoupled mcMove = new MCMoveAtomCoupled(potentialMaster, meterPE, getRandom(), space);
        mcMove.setPotential(p2);
        integrator.getMoveManager().addMCMove(mcMove);

        box.getBoundary().getEventManager().removeListener(potentialMaster.getNbrCellManager(box));

        if (!nvt) {
            // using Carnahan-Starling EOS.  the pressure will be too high because
            // we have a solid, but OK.
            double eta = rho * Math.PI / 6;
            double den = 1 - eta;
            double z = (1 + eta + eta * eta - eta * eta * eta) / (den * den * den);
            double p = z * rho;

            // hard coded pressure for rho=1.2
            if (space.D() == 2) {
                p = 14.9;
            } else {
                if (rho == 1.2) {
                    p = 23.3;
                } else {
                    p = 47.6;
                }
            }

            MCMove mcMoveVolume;
            if (fancyMove) {
                // fancy move
                mcMoveVolume = new MCMoveVolumeSolid(potentialMaster, coordinateDefinition, getRandom(), space, p);
                ((MCMoveVolumeSolid) mcMoveVolume).setTemperature(1.0);
            } else {
                // standard move
                mcMoveVolume = new MCMoveVolume(potentialMaster, getRandom(), space, p);
            }
            ((MCMoveStepTracker) mcMoveVolume.getTracker()).setNoisyAdjustment(true);
            integrator.getMoveManager().addMCMove(mcMoveVolume);
        }
    }
    
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        final HSMD3DParameters params = new HSMD3DParameters();
        ParseArgs parseArgs = new ParseArgs(params);
        parseArgs.parseArgs(args);
        final HSNPT sim = new HSNPT(params.numAtoms == 120 ? Space2D.getInstance() : Space3D.getInstance(), params.numAtoms, params.rho, params.nvt, params.fancyMove, params.sigma2);
        
        DataFork volumeFork = null;
        AccumulatorAverageFixed volumeAvg = null;
        final FileWriter vfw;
        AccumulatorAverageFixed displacementAvg = null;
        AccumulatorAverageFixed displacementMax = null;
        AccumulatorAverageFixed maxExpansionAvg = null;
        if (!params.nvt) {
            MeterVolume meterVolume = new MeterVolume();
            meterVolume.setBox(sim.box);
            volumeFork = new DataFork();
            DataPumpListener volumePump = new DataPumpListener(meterVolume, volumeFork, params.numAtoms);
            sim.integrator.getEventManager().addListener(volumePump);
            long nData = params.numSteps/params.numAtoms;
            long blockSize = (nData+99)/100;
            volumeAvg = new AccumulatorAverageFixed(blockSize);
            volumeFork.addDataSink(volumeAvg);
        }
        else {
//            MeterDisplacement meterDisplacement = new MeterDisplacement(sim.space, sim.coordinateDefinition, 0.001);
//            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(meterDisplacement, params.numAtoms));
    
            MeterDisplacementRMS meterDisplacementRMS = new MeterDisplacementRMS(sim.space, sim.coordinateDefinition);
            displacementAvg = new AccumulatorAverageFixed();
            DataPumpListener displacementAvgPump = new DataPumpListener(meterDisplacementRMS, displacementAvg, params.numAtoms);
            sim.integrator.getEventManager().addListener(displacementAvgPump);
    
//            MeterDisplacementMax meterDisplacementMax = new MeterDisplacementMax(sim.space, sim.coordinateDefinition, 0.001);
//            displacementMax = new AccumulatorAverageFixed();
//            DataPumpListener displacementMaxPump = new DataPumpListener(meterDisplacementMax, displacementMax, params.numAtoms);
//            sim.integrator.getEventManager().addListener(displacementMaxPump);
    
            MeterMaxExpansion meterMaxExpansion = new MeterMaxExpansion(sim.space, sim.box, sim.potentialMaster.getNeighborManager(sim.box));
            maxExpansionAvg = new AccumulatorAverageFixed();
            DataPumpListener maxExpansionPump = new DataPumpListener(meterMaxExpansion, maxExpansionAvg, params.numAtoms);
            sim.integrator.getEventManager().addListener(maxExpansionPump);
        }

        if (true) {
            SimulationGraphic graphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE);

            if (!params.nvt) {
                AccumulatorHistory densityHistory = new AccumulatorHistory(new HistoryCollapsingAverage());
//                densityHistory.setPushInterval();
                densityHistory.setTimeDataSource(new DataSourceCountSteps(sim.integrator));
                volumeFork.addDataSink(densityHistory);
                DisplayPlot densityPlot = new DisplayPlot();
                densityHistory.setDataSink(densityPlot.getDataSet().makeDataSink());
                densityPlot.setLabel("volume");
                graphic.add(densityPlot);
            }
            
            ActionSummer summer = new ActionSummer(sim.box, sim.getSpace());
            sim.integrator.getEventManager().addListener(new IntegratorListenerAction(summer, params.numAtoms));
            final ColorSchemeDisplacement colorScheme = new ColorSchemeDisplacement(sim.coordinateDefinition, sim.getSpace(), summer);
            graphic.getDisplayBox(sim.box).setColorScheme(colorScheme);
            
            final DiameterHashByType diameterHash = (DiameterHashByType)graphic.getDisplayBox(sim.box).getDiameterHash();
            diameterHash.setDiameter(sim.getSpecies(0).getAtomType(0), 1.0);
            if (params.sigma2 != 1) {
                diameterHash.setDiameter(sim.getSpecies(1).getAtomType(0), params.sigma2);
            }
            
            DeviceSlider siteFactorSlider = new DeviceSlider(sim.getController(), new Modifier() {
                
                public void setValue(double newValue) {
                    if (newValue < 0 || newValue > 6) throw new IllegalArgumentException();
                    double x = Math.pow(10, newValue);
                    colorScheme.setSiteFactor(x);
                }
                
                public double getValue() {
                    return Math.log10(colorScheme.getSiteFactor());
                }
                
                public String getLabel() {
                    return "site factor";
                }
                
                public Dimension getDimension() {
                    return Null.DIMENSION;
                }
            });
            siteFactorSlider.setPrecision(1);
            siteFactorSlider.setMinimum(1);
            siteFactorSlider.setMaximum(4);
            siteFactorSlider.setLabel("Site factor");
            siteFactorSlider.setShowBorder(true);
            graphic.add(siteFactorSlider);
            
            DeviceSlider fluctuationFactorSlider = new DeviceSlider(sim.getController(), new Modifier() {
                
                public void setValue(double newValue) {
                    if (newValue < 0 || newValue > 6) throw new IllegalArgumentException();
                    double x = Math.pow(10, newValue);
                    colorScheme.setFluctuationFactor(x);
                }
                
                public double getValue() {
                    return Math.log10(colorScheme.getFluctuationFactor());
                }
                
                public String getLabel() {
                    return "fluctuation factor";
                }
                
                public Dimension getDimension() {
                    return Null.DIMENSION;
                }
            });
            fluctuationFactorSlider.setPrecision(1);
            fluctuationFactorSlider.setMinimum(1);
            fluctuationFactorSlider.setMaximum(4);
            fluctuationFactorSlider.setLabel("Fluctuation factor");
            fluctuationFactorSlider.setShowBorder(true);
            graphic.add(fluctuationFactorSlider);
            
            DeviceBox nominalFluctuationBox = new DeviceBox();
            nominalFluctuationBox.setController(sim.getController());
            nominalFluctuationBox.setEditable(true);
            nominalFluctuationBox.setLabel("Nominal fluctuation");
            nominalFluctuationBox.setModifier(new Modifier() {
                
                public void setValue(double newValue) {
                    if (newValue < 0 || newValue > 1) throw new IllegalArgumentException();
                    colorScheme.setNominalFluctuation(newValue);
                }
                
                public double getValue() {
                    return colorScheme.getNominalFluctuation();
                }
                
                public String getLabel() {
                    return "Nominal fluctuation";
                }
                
                public Dimension getDimension() {
                    return Length.DIMENSION;
                }
            });
            graphic.add(nominalFluctuationBox);
            
            
            DeviceSelector fluctuationSelector = new DeviceSelector(sim.getController());
            fluctuationSelector.setLabel("Fluctuation about");
            fluctuationSelector.addOption("lattice site", new IAction() {
                public void actionPerformed() {
                    colorScheme.setIsFluctuationFromAvg(false);
                }
            });
            fluctuationSelector.addOption("avg position", new IAction() {
                public void actionPerformed() {
                    colorScheme.setIsFluctuationFromAvg(true);
                }
            });
            fluctuationSelector.setSelected(1);

            graphic.add(fluctuationSelector);
            
            DeviceBox sigmaBox = new DeviceBox();
            sigmaBox.setController(sim.getController());
            sigmaBox.setModifier(new Modifier() {

                public double getValue() {
                    return (sim.pCross.getCollisionDiameter() * 2 - 1);
                }
                
                public void setValue(double newValue) {
                    double oldValue = getValue();
                    if  (newValue < 0.5 || newValue > 2) throw new IllegalArgumentException();
                    sim.pCross.setCollisionDiameter((newValue+1)*0.5);
                    try {
                        sim.integrator.reset();
                    }
                    catch (ConfigurationOverlapException e) {
                        newValue = oldValue;
                        sim.pCross.setCollisionDiameter((newValue+1)*0.5);
                        sim.integrator.reset();
                        return;
                    }
                    diameterHash.setDiameter(sim.getSpecies(1).getAtomType(0), newValue);
                }
                
                public String getLabel() {
                    return "Diameter";
                }
                
                public Dimension getDimension() {
                    return Length.DIMENSION;
                }
            });
            sigmaBox.setLabel("Sigma");
            graphic.add(sigmaBox);
            
            graphic.makeAndDisplayFrame();
            
            return;
        }
        long t1 = System.currentTimeMillis();
        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numSteps/10);
        if (!params.nvt) {
            volumeAvg.reset();
        }
        else {
            displacementAvg.reset();
            displacementMax.reset();
            maxExpansionAvg.reset();
        }
            
        System.out.println("equilibration finished");
        sim.integrator.getMoveManager().setEquilibrating(false);

        if (!params.nvt && false) {
            final String vfile = "npt_"+(params.fancyMove ? "fancy" : "simple")+"_vt.dat";
            try {
                vfw = new FileWriter(vfile);
            }
            catch (IOException e) {
                throw new RuntimeException(e);
            }
            volumeFork.addDataSink(new IDataSink() {
                
                public void putDataInfo(IDataInfo dataInfo) {}
                
                public void putData(IData data) {
                    try {
                        long iStep = sim.integrator.getStepCount()/params.numAtoms;
                        if (iStep % 10 == 0 || params.fancyMove) {
                            vfw.write(iStep +" "+data.getValue(0)+"\n");
                        }
                    }
                    catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }
            });
        }
        else {
            vfw = null;
        }

        sim.getController().runActivityBlocking(new ActivityIntegrate(sim.integrator), params.numSteps);
        System.out.println("time "+(System.currentTimeMillis()-t1)/1000);

        if (!params.nvt) {
            if (vfw != null) {
                try {
                    vfw.close();
                }
                catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            double vavg = volumeAvg.getData().getValue(volumeAvg.AVERAGE.index);
            double verr = volumeAvg.getData().getValue(volumeAvg.ERROR.index);
            double vstdev = volumeAvg.getData().getValue(volumeAvg.STANDARD_DEVIATION.index);
            double vcorr = volumeAvg.getData().getValue(volumeAvg.BLOCK_CORRELATION.index);
            System.out.println("avg volume "+vavg+"  err "+verr+"  stdev "+vstdev+"  correlation "+vcorr);
            System.out.println("avg density "+params.numAtoms/vavg+" "+params.numAtoms/(vavg*vavg)*verr);
        }

        if (params.nvt) {
            double davg = displacementAvg.getData().getValue(displacementAvg.AVERAGE.index);
            double dstdev = displacementAvg.getData().getValue(displacementAvg.STANDARD_DEVIATION.index);
            double derr = displacementAvg.getData().getValue(displacementAvg.ERROR.index);
            System.out.println("displacement avg "+davg+" stdev "+dstdev+" err "+derr);

//            double dmaxavg = displacementMax.getData().getValue(displacementMax.StatType.AVERAGE.index);
//            double dmaxstdev = displacementMax.getData().getValue(displacementMax.StatType.STANDARD_DEVIATION.index);
//            System.out.println("displacement max avg "+dmaxavg+" stdev "+dmaxstdev);

            double emaxavg = maxExpansionAvg.getData().getValue(maxExpansionAvg.AVERAGE.index);
            double emaxstdev = maxExpansionAvg.getData().getValue(maxExpansionAvg.STANDARD_DEVIATION.index);
            System.out.println("max expansion avg "+emaxavg+" stdev "+emaxstdev);
        }
    }
    
    public static class ColorSchemeDisplacement extends ColorScheme {
        protected final Vector v, v2;
        protected final ActionSummer summer;
        protected CoordinateDefinition coordinateDefinition;
        protected AtomLeafAgentManager agentManager;
        protected double siteFactor = 1000;
        protected double fluctuationFactor = 1000;
        protected double nominalFluctuation = 0.0627;  // 256=>0.0607, 864=>0.0627
        protected boolean fluctuationFromAvg = true;
        protected double mysum;
        
        public ColorSchemeDisplacement(CoordinateDefinition cDef, Space space, ActionSummer summer) {
            super();
            v = space.makeVector();
            v2 = space.makeVector();
            coordinateDefinition = cDef;
            agentManager = summer.getAgentManager();
            this.summer = summer;
        }
        
        public double getSiteFactor() {
            return siteFactor;
        }
        
        public void setSiteFactor(double newSiteFactor) {
            siteFactor = newSiteFactor;
        }
        
        public double getFluctuationFactor() {
            return fluctuationFactor;
        }
        
        public void setFluctuationFactor(double newFluctuationFactor) {
            fluctuationFactor = newFluctuationFactor;
        }
        
        public boolean getIsFluctuationFromAvg() {
            return fluctuationFromAvg;
        }
        
        public void setIsFluctuationFromAvg(boolean newFluctuationFromAvg) {
            fluctuationFromAvg = newFluctuationFromAvg;
        }
        
        public double getNominalFluctuation() {
            return nominalFluctuation;
        }
        
        public void setNominalFluctuation(double newNominalFlucutation) {
            nominalFluctuation = newNominalFlucutation;
        }
        
        public Color getAtomColor(IAtom a) {
            MyAgent agent = (MyAgent)agentManager.getAgent(a);
            int count = summer.getCount();
            v.Ea1Tv1(1.0/count, agent.pSum);
            int red = (int)Math.round(255*(siteFactor * (v.Mv1Squared(coordinateDefinition.getLatticePosition(a)))));
            red = (red / 10) * 10;
            if (red > 255) red = 255;
            if (fluctuationFromAvg) {
                v.TE(v);
                v.TE(-1);
                v.PEa1Tv1(1.0/count, agent.pSum2);
            }
            else {
                v.E(coordinateDefinition.getLatticePosition(a));
                v.PEa1Tv1(-2.0/count, agent.pSum);
                v.TE(coordinateDefinition.getLatticePosition(a));
                v.PEa1Tv1(1.0/count, agent.pSum2);
            }
            double fsum = Math.sqrt(v.getX(0) + v.getX(1) + v.getX(2)) - nominalFluctuation;
//            mysum += Math.sqrt(v.getX(0) + v.getX(1) + v.getX(2));
//            if (a.getLeafIndex() == 255) {
//                System.out.println(mysum/256);
//                mysum = 0;
//            }
            int c = (int)Math.round(255*(fluctuationFactor * fsum));
//            if (a.getLeafIndex()==0) System.out.println(fsum);
            c = (c / 10) * 10;
            int blue = 0, green = 0;
            if (c>0) {
                blue = c;
                if (blue > 255) blue = 250;
            }
            else {
                green = -c;
                if (green > 255) green = 250;
            }
            return new Color(red, green, blue);
        }
    }
    
    public static class ActionSummer implements IAction, AgentSource<MyAgent> {
        
        protected final AtomLeafAgentManager<MyAgent> agentManager;
        protected final Vector v;
        protected final Space space;
        protected int count;
        
        public ActionSummer(Box box, Space space) {
            this.space = space;
            agentManager = new AtomLeafAgentManager<MyAgent>(this, box);
            v = space.makeVector();
        }
        
        public AtomLeafAgentManager<MyAgent> getAgentManager() {
            return agentManager;
        }
        
        public void actionPerformed() {
            Box box = agentManager.getBox();
            IAtomList atoms = box.getLeafList();
            for (int i = 0; i<atoms.size(); i++) {
                IAtom atom = atoms.get(i);
                MyAgent agent = agentManager.getAgent(atom);
                v.E(atom.getPosition());
                agent.pSum.PE(v);
                v.TE(v);
                agent.pSum2.PE(v);
            }
            count++;
        }
        
        public int getCount() {
            return count;
        }
        
        public void reset() {
            Box box = agentManager.getBox();
            IAtomList atoms = box.getLeafList();
            for (int i = 0; i<atoms.size(); i++) {
                IAtom atom = atoms.get(i);
                MyAgent agent = agentManager.getAgent(atom);
                agent.pSum.E(0);
                agent.pSum2.E(0);
            }
            count = 0;
        }
        
        public MyAgent makeAgent(IAtom a, Box agentBox) {
            return new MyAgent(space);
        }
        
        public void releaseAgent(MyAgent agent, IAtom atom, Box agentBox) {}

    }
    
    public static class MyAgent {
        public Vector pSum, pSum2;
        public MyAgent(Space space) {
            pSum = space.makeVector();
            pSum2 = space.makeVector();
        }
    }
    
    public static class HSMD3DParameters extends ParameterBase {
        public int numAtoms = 864;
        public double rho = 1.2;
        public long numSteps = 10000000;
        public boolean nvt = true;
        public boolean fancyMove = false;
        public double sigma2 = 1.11;//1.11;
    }
}
