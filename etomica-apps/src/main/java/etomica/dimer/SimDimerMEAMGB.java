/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.action.CalcVibrationalModes;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.chem.elements.Tin;
import etomica.config.ConfigurationFile;
import etomica.config.GrainBoundaryTiltConfiguration;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.history.HistoryCollapsingAverage;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.numerical.CalcGradientDifferentiable;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionTypesCombination;
import etomica.nbr.list.PotentialMasterList;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class SimDimerMEAMGB extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerMEAMadatomSn";
    public final PotentialMasterList potentialMaster;
    public final PotentialMasterListDimer potentialMasterD;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorDimerMin;
    public Box box;
    public double clusterRadius;
    public Vector[] saddle;
    public SpeciesSpheresMono fixed, movable, dimer;
    public PotentialMEAM potential;
    public PotentialCalculationForcePressureSumGB pcGB;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    public CalcGradientDifferentiable calcGradientDifferentiable;
    public CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    public Vector adAtomPos;
    public IMoleculeList movableSet;
    public int [] millerPlane;
    

    
    public SimDimerMEAMGB(int[] amillerPlane, int[] boxSize) {
    	super(Space3D.getInstance());
    	
    	this.millerPlane = amillerPlane;
    	potentialMaster = new PotentialMasterList(this, space);
    	potentialMasterD = new PotentialMasterListDimer(this, space);
        
      //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(2, 5, space), space);
        addBox(box);
     
      //SPECIES
        
        //Sn
        Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);
        Tin dimerTin = new Tin("SnD", 118.710);
        fixed = new SpeciesSpheresMono(space, tinFixed);
        fixed.setIsDynamic(true);
        movable = new SpeciesSpheresMono(space, Tin.INSTANCE);
        fixed.setIsDynamic(true);
        dimer = new SpeciesSpheresMono(space, dimerTin);
        fixed.setIsDynamic(true);
        addSpecies(fixed);
        addSpecies(movable);
        addSpecies(dimer);
        potential = new PotentialMEAM(space);
        potential.setParameters(fixed.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(movable.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(dimer.getLeafType(), ParameterSetMEAM.Sn);
        
        
        
        //Sn
        //beta-Sn box
        
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
              
        double a = 5.92; 
        double c = 3.23;
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, a, c);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new ISpecies[] {fixed, movable}, potential.getRange(), space);
            
        
        //Ag
        /**
        Silver silverFixed = new Silver("AgFix", Double.POSITIVE_INFINITY);
        fixed = new SpeciesSpheresMono(this, Silver.INSTANCE);
        movable = new SpeciesSpheresMono(this, Silver.INSTANCE);
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        ((AtomTypeSphere)fixed.getLeafType()).setDiameter(2.8895); 
        ((AtomTypeSphere)movable.getLeafType()).setDiameter(2.8895);
        potential = new PotentialMEAM(space);
        potential.setParameters(agFix, ParameterSetMEAM.Ag);
        potential.setParameters(ag, ParameterSetMEAM.Ag);
        potential.setParameters(agAdatom, ParameterSetMEAM.Ag);
        potential.setParameters(movable, ParameterSetMEAM.Ag);
        
        double a = 4.0863;
        PrimitiveCubic primitive = new PrimitiveCubic(space, a);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new ISpecies[] {fixed, movable}, 4.56, space);

        */
    
        
        //Cu
       /**
        //Copper copperFixed = new Copper("CuFix", Double.POSITIVE_INFINITY);
        fixed = new SpeciesSpheresMono(this, space, Copper.INSTANCE);
        movable = new SpeciesSpheresMono(this, space, Copper.INSTANCE);
        dimer = new SpeciesSpheresMono(this, space, Copper.INSTANCE);
        getSpeciesManager().addSpecies(fixed);
        getSpeciesManager().addSpecies(movable);
        getSpeciesManager().addSpecies(dimer);
        ((AtomTypeSphere)fixed.getLeafType()).setDiameter(2.5561); 
        ((AtomTypeSphere)dimer.getLeafType()).setDiameter(2.5561); 
        ((AtomTypeSphere)movable.getLeafType()).setDiameter(2.5561);
        potential = new PotentialMEAM(space);
        potential.setParameters(fixed.getLeafType(), ParameterSetMEAM.Cu);
        potential.setParameters(movable.getLeafType(), ParameterSetMEAM.Cu);
        potential.setParameters(dimer.getLeafType(), ParameterSetMEAM.Cu);
        
        double a = 3.6148;
        PrimitiveCubic primitive = new PrimitiveCubic(space, a);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());
        GrainBoundaryTiltConfiguration gbtilt = new GrainBoundaryTiltConfiguration(crystal, crystal, new ISpecies[] {fixed, movable}, potential.getRange(), space);
       */

        this.potentialMaster.addPotential(potential, new AtomType[]{fixed.getLeafType(), movable.getLeafType(), dimer.getLeafType()});
        potentialMaster.setRange(potential.getRange()*1.1);
        CriterionSimple criteria = new CriterionSimple(this, space, potential.getRange(), potential.getRange()*1.1);
        potentialMaster.setCriterion(potential, new CriterionTypesCombination(criteria, new AtomType[]{fixed.getLeafType(), movable.getLeafType(), dimer.getLeafType()}));

        this.potentialMasterD.addPotential(potential, new AtomType[]{movable.getLeafType(), dimer.getLeafType()});
        potentialMasterD.setSpecies(new ISpecies []{dimer, movable});
        potentialMasterD.setRange(potential.getRange()*1.1);
        CriterionSimple criteria2 = new CriterionSimple(this, space, potential.getRange(), potential.getRange()*1.1);
        potentialMasterD.setCriterion(potential, new CriterionTypesCombination(criteria2, new AtomType[]{movable.getLeafType(), dimer.getLeafType()}));
        
        gbtilt.setFixedSpecies(fixed);
        gbtilt.setMobileSpecies(movable);

        gbtilt.setGBplane(millerPlane);
        gbtilt.setBoxSize(box, boxSize);
        gbtilt.initializeCoordinates(box);
               
        Vector newBoxLength = space.makeVector();
        newBoxLength.E(box.getBoundary().getBoxSize());
        newBoxLength.setX(2,newBoxLength.getX(2)+1.0);
        newBoxLength.setX(1,newBoxLength.getX(1)+0.0001);
        newBoxLength.setX(0,newBoxLength.getX(0)+0.0001);
        box.getBoundary().setBoxSize(newBoxLength);
        
        
    }

    public static void main(String[] args) {
        String fileName = "sngb101v-test";//args[0];
        //int mdSteps = 10;//Integer.parseInt(args[1]);
        /*
        int h = Integer.parseInt(args[1]);
        int k = Integer.parseInt(args[2]);
        int l = Integer.parseInt(args[3]);

        int x = Integer.parseInt(args[4]);
        int y = Integer.parseInt(args[5]);
        int z = Integer.parseInt(args[6]);
        */
        final String APP_NAME = "SimDimerMEAMGBCluster";

        final SimDimerMEAMGB sim = new SimDimerMEAMGB(new int[]{2, 1, 0}, new int[]{2, 6, 12});

        sim.initializeConfiguration("sngb210-2612");

        //System.out.println(sim.box.getBoundary().getDimensions().get(0));


        Vector dimerCenter = sim.getSpace().makeVector();
        dimerCenter.setX(0, sim.box.getBoundary().getBoxSize().getX(0) / 2.0);
        dimerCenter.setX(1, 1.0);
        dimerCenter.setX(2, 0.0);
        Vector cubeSize = sim.getSpace().makeVector();
        cubeSize.setX(0, 6.0);
        cubeSize.setX(1, 8.0);
        cubeSize.setX(2, 8.0);

        if (sim.millerPlane[2] == 0) {
            dimerCenter.setX(1, sim.box.getBoundary().getBoxSize().getX(1) / 2.0);
            dimerCenter.setX(0, 1.0);
            dimerCenter.setX(2, 0.0);
            cubeSize.setX(0, 6.0);
            cubeSize.setX(1, 8.0);
            cubeSize.setX(2, 8.0);
        }

        IAtomList list = sim.box.getLeafList();
        Vector rij = sim.space.makeVector();
        Vector3D move = new Vector3D(0.0, 0.0, 5.0);
        Vector3D move2 = new Vector3D(0.0, 0.0, 10.0);
        move2.PE(sim.box.getBoundary().getBoxSize());
        System.out.println("Atoms: " + list.getAtomCount());
        System.out.println("Interface Area: " + move2.getX(0) * move2.getX(1) + " angstroms");

        /*
        sim.box.getBoundary().setDimensions(move2);
        for(int i=0; i<list.getAtomCount(); i++){
        	rij = list.getAtom(i)).getPosition();
        	if(rij.get(2)>0.0){rij.PE(move);}
        	else{rij.ME(move);}
        }
        */
        //sim.setMovableAtomsSphere(6.0, dimerCenter);
        sim.setMovableAtomsCube(cubeSize, dimerCenter);
        sim.setMovableAtomsList();

        sim.removeAtoms(3.0, dimerCenter);
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */

        //sim.initializeConfiguration("sngb101-d1-00_saddle");
        sim.enableMolecularDynamics(2);
        //sim.enableDimerSearch(fileName, 2000, false, false);
        //sim.integratorDimer.setRotNum(1);
        //sim.enableMinimumSearch("sngb5-r1", true);

        /*
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"_saddle.getyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimer.addIntervalAction(xyzwriter);
        sim.integratorDimer.setActionInterval(xyzwriter, 5);
        */


        /*
        WriteConfiguration writer = new WriteConfiguration(sim.getSpace());
        writer.setBox(sim.box);
        writer.setConfName(fileName);
        sim.integratorMD.addIntervalAction(writer);
        sim.integratorMD.setActionInterval(writer, 10000);
        */

        MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMasterD);
        energyMeter.setBox(sim.box);
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorAverageCollapsing accumulatorAveragePE = new AccumulatorAverageCollapsing();
        DataPump energyPump = new DataPump(energyMeter, accumulatorAveragePE);
        accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{AccumulatorAverage.MOST_RECENT});
        DisplayPlot plotPE = new DisplayPlot();
        plotPE.setLabel("PE Plot");
        energyAccumulator.setDataSink(plotPE.getDataSet().makeDataSink());
        accumulatorAveragePE.setPushInterval(1);
        //sim.integratorDimer.addIntervalAction(energyPump);
        //sim.integratorDimer.setActionInterval(energyPump,1);


        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1);
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(plotPE);

        sim.integratorMD.getEventManager().addListener(new IntegratorListenerAction(simGraphic.getPaintAction(sim.box)));

        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());

        colorScheme.setColor(sim.fixed.getLeafType(), java.awt.Color.blue);
        colorScheme.setColor(sim.movable.getLeafType(), java.awt.Color.gray);
        colorScheme.setColor(sim.dimer.getLeafType(), java.awt.Color.orange);

        simGraphic.makeAndDisplayFrame(APP_NAME);
    }
    
    public void setMovableAtomsSphere(double distance, Vector center){
        distance = distance*distance;
        Vector rij = space.makeVector();
        MoleculeArrayList movableList = new MoleculeArrayList();
        IMoleculeList loopSet = box.getMoleculeList();
        for (int i=0; i<loopSet.getMoleculeCount(); i++){
            if(loopSet.getMolecule(i).getType()==fixed){continue;}
        	rij.E(loopSet.getMolecule(i).getChildList().getAtom(0).getPosition());
            rij.Ev1Mv2(center, rij);
            box.getBoundary().nearestImage(rij);
            if(rij.squared()<distance){//Math.abs(rij.get(0)) < 0.5 && Math.abs(rij.get(1)) < distance && Math.abs(rij.get(2)) < distance){
               movableList.add(loopSet.getMolecule(i));
            }
        }
        for (int i=0; i<movableList.getMoleculeCount(); i++){
            IMolecule newMolecule = dimer.makeMolecule();
            box.addMolecule(newMolecule);
           newMolecule.getChildList().getAtom(0).getPosition().E(movableList.getMolecule(i).getChildList().getAtom(0).getPosition());
           box.removeMolecule(movableList.getMolecule(i));
        }
        movableSet = box.getMoleculeList(dimer);
    }
    
    public void setMovableAtomsCube(Vector dimensions, Vector center){
        Vector cube = dimensions;
        Vector rij = space.makeVector();
        MoleculeArrayList movableList = new MoleculeArrayList();
        IMoleculeList loopSet = box.getMoleculeList();
        for (int i=0; i<loopSet.getMoleculeCount(); i++){
        	if(loopSet.getMolecule(i).getType()==fixed){continue;}
            rij.E(loopSet.getMolecule(i).getChildList().getAtom(0).getPosition());
            rij.Ev1Mv2(center, rij);
            box.getBoundary().nearestImage(rij);
            if(Math.abs(rij.getX(0)) < cube.getX(0) && Math.abs(rij.getX(1)) < cube.getX(1) && Math.abs(rij.getX(2)) < cube.getX(2)){
               movableList.add(loopSet.getMolecule(i));
            }
        }
        for (int i=0; i<movableList.getMoleculeCount(); i++){
            IMolecule newMolecule = dimer.makeMolecule();
            box.addMolecule(newMolecule);
           newMolecule.getChildList().getAtom(0).getPosition().E(movableList.getMolecule(i).getChildList().getAtom(0).getPosition());
           box.removeMolecule(movableList.getMolecule(i));
        }
        movableSet = box.getMoleculeList(dimer);
    }
    
    public void setMovableAtomsList(){
        MoleculeArrayList neighborList = new MoleculeArrayList();
        MoleculeArrayList fixedList = new MoleculeArrayList();
        IMoleculeList loopSet = box.getMoleculeList();
        IMoleculeList dimerSet = box.getMoleculeList(dimer);
        for(int i=0; i<loopSet.getMoleculeCount(); i++){
            if(loopSet.getMolecule(i).getType()==dimer){
                continue;
            }
            boolean fixedFlag = true;
            for(int j=0; j<dimerSet.getMoleculeCount(); j++){
                Vector dist = space.makeVector();
                dist.Ev1Mv2(loopSet.getMolecule(i).getChildList().getAtom(0).getPosition(),dimerSet.getMolecule(j).getChildList().getAtom(0).getPosition());
                box.getBoundary().nearestImage(dist);
                if(Math.sqrt(dist.squared())<potentialMasterD.getMaxPotentialRange()+2.0){
                    neighborList.add(loopSet.getMolecule(i));
                    fixedFlag = false;
                    break;
                }

            }
            if(fixedFlag){
                fixedList.add(loopSet.getMolecule(i));
            }
        }
        for (int i=0; i<neighborList.getMoleculeCount(); i++){
            IMolecule newMolecule = movable.makeMolecule();
            box.addMolecule(newMolecule);
            newMolecule.getChildList().getAtom(0).getPosition().E(neighborList.getMolecule(i).getChildList().getAtom(0).getPosition());
            box.removeMolecule(neighborList.getMolecule(i));
         }
        for (int i=0; i<fixedList.getMoleculeCount(); i++){
            IMolecule newMolecule = fixed.makeMolecule();
            box.addMolecule(newMolecule);
            newMolecule.getChildList().getAtom(0).getPosition().E(fixedList.getMolecule(i).getChildList().getAtom(0).getPosition());
            box.removeMolecule(fixedList.getMolecule(i));
         }

    }
    
    //Must be run after setMovableAtoms
    public void removeAtoms(double distance, Vector center){
        distance = distance*distance;
        int rmvCount = 0;
        Vector rij = space.makeVector();
        //movable species
        IMoleculeList loopSet = box.getMoleculeList(movable);
        for (int i=0; i<loopSet.getMoleculeCount(); i++){
            rij.Ev1Mv2(center,loopSet.getMolecule(i).getChildList().getAtom(0).getPosition());
            box.getBoundary().nearestImage(rij);
            if(rij.squared() < distance){
               box.removeMolecule(loopSet.getMolecule(i));
               rmvCount++;
            }
        }

        //dimer species
        IMoleculeList loopSet2 = box.getMoleculeList(dimer);
        for (int i=0; i<loopSet2.getMoleculeCount(); i++){
            rij.Ev1Mv2(center,loopSet2.getMolecule(i).getChildList().getAtom(0).getPosition());
            box.getBoundary().nearestImage(rij);
            if(rij.squared() < distance){
               box.removeMolecule(loopSet2.getMolecule(i));
               rmvCount++;
            }
        }
        System.out.println(rmvCount+" atoms removed.");
    }
    
    public void initializeConfiguration(String fileName){
        ConfigurationFile config = new ConfigurationFile(fileName);
        config.initializeCoordinates(box);
    }

    public void generateConfigs(String fileName, double percentd) {

        Vector workVector = space.makeVector();
        Vector[] currentPos = new Vector[movableSet.getMoleculeCount()];
        for(int i=0; i<currentPos.length; i++){
            currentPos[i] = space.makeVector();
            currentPos[i].E(movableSet.getMolecule(i).getChildList().getAtom(0).getPosition());
        }

        //Create multiple configurations
        for(int m=0; m<50; m++){
            WriteConfiguration genConfig = new WriteConfiguration(space);
            genConfig.setBox(box);
            genConfig.setConfName(fileName+"_config_"+m);
            //Displaces atom's by at most +/-0.03 in each coordinate
            for(int i=0; i<movableSet.getMoleculeCount(); i++){
                Vector atomPosition = movableSet.getMolecule(i).getChildList().getAtom(0).getPosition();
                for(int j=0; j<3; j++){
                    workVector.setX(j,percentd*random.nextGaussian());
                }
                atomPosition.Ev1Pv2(currentPos[i],workVector);
            }
            genConfig.actionPerformed();
        }
    }

    public void enableMolecularDynamics(long maxSteps) {
        integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, space, box);
        integratorMD.setTimeStep(0.001);
        integratorMD.setTemperature(Kelvin.UNIT.toSim(100));
        integratorMD.setThermostatInterval(100);
        integratorMD.setIsothermal(true);
        //pcGB = new PotentialCalculationForcePressureSumGB(space, box);
        //integratorMD.setForceSum(pcGB);
        integratorMD.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        activityIntegrateMD = new ActivityIntegrate(integratorMD);
        getController().addAction(activityIntegrateMD);
        activityIntegrateMD.setMaxSteps(maxSteps);
    }

    public void enableDimerSearch(String fileName, long maxSteps, Boolean orthoSearch, Boolean fine) {

        integratorDimer = new IntegratorDimerRT(this, potentialMasterD, new ISpecies[]{dimer}, space, box);
        integratorDimer.setOrtho(orthoSearch, false);
        if (fine) {
            ConfigurationFile configFile = new ConfigurationFile(fileName + "_saddle");
            configFile.initializeCoordinates(box);

            integratorDimer.setFileName(fileName + "_fine");
            integratorDimer.deltaR = 0.0005;
            integratorDimer.dXl = 10E-5;
            integratorDimer.deltaXmax = 0.005;
            integratorDimer.dFsq = 0.0001 * 0.0001;
            integratorDimer.dFrot = 0.01;
        }
        integratorDimer.setFileName(fileName);
        integratorDimer.getEventManager().addListener(potentialMasterD.getNeighborManager(box));
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
        integratorDimer.setActivityIntegrate(activityIntegrateDimer);
        getController().addAction(activityIntegrateDimer);
        activityIntegrateDimer.setMaxSteps(maxSteps);
    }

    public void enableMinimumSearch(String fileName, Boolean normalDir) {

        integratorDimerMin = new IntegratorDimerMin(this, potentialMasterD, new ISpecies[]{dimer}, normalDir, space, box);
        integratorDimerMin.getEventManager().addListener(potentialMasterD.getNeighborManager(box));
        activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
        integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
        getController().addAction(activityIntegrateMin);
    }
    
}
