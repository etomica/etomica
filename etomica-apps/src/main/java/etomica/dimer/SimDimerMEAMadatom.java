/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.action.CalcVibrationalModes;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.atom.AtomType;
import etomica.atom.MoleculeArrayList;
import etomica.box.Box;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
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
import etomica.listener.IntegratorListenerAction;
import etomica.math.numerical.CalcGradientDifferentiable;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionTypesCombination;
import etomica.nbr.list.PotentialMasterList;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class SimDimerMEAMadatom extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerMEAMadatomSn";
    public final PotentialMasterList potentialMaster;
    public final PotentialMasterListDimer potentialMasterD;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorDimerMin;
    public Box box;
    public Vector[] saddle, normal;
    public SpeciesSpheresMono fixed, potentialSpecies, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    public CalcGradientDifferentiable calcGradientDifferentiable;
    public CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    public Vector adAtomPos;
    public IMoleculeList movableSet;
    //public Boolean saddleFine, calcModes, minSearch, normalDir;
    
    public SimDimerMEAMadatom() {
        super(Space3D.getInstance());    	
        potentialMaster = new PotentialMasterList(this, space);
        potentialMasterD = new PotentialMasterListDimer(this, space);
        
      //SIMULATION BOX
        box = new Box(new BoundaryRectangularSlit(0, 5, space), space);
        addBox(box);
     
      //SPECIES
        
        //Sn
        Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);  
        fixed = new SpeciesSpheresMono(space, tinFixed);
        fixed.setIsDynamic(true);
        movable = new SpeciesSpheresMono(space, Tin.INSTANCE);
        movable.setIsDynamic(true);
        potentialSpecies = new SpeciesSpheresMono(space, tinFixed);
        potentialSpecies.setIsDynamic(true);
        addSpecies(fixed);
        addSpecies(movable);
        addSpecies(potentialSpecies);
        box.setNMolecules(fixed, 180);
        
        potential = new PotentialMEAM(space);
        potential.setParameters(fixed.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(movable.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(potentialSpecies.getLeafType(), ParameterSetMEAM.Sn);
        
        
        
        //Sn
        //beta-Sn box
        
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta box (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
              
        double a = 5.92; 
        double c = 3.23;
        box.getBoundary().setBoxSize(new Vector3D(a*3, a*3, c*5));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, a, c);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());

                
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
        */
        
        Configuration config = new ConfigurationLattice(crystal, space);
        config.initializeCoordinates(box);

        this.potentialMaster.addPotential(potential, new AtomType[]{movable.getLeafType(), potentialSpecies.getLeafType()});
        potentialMaster.setRange(potential.getRange()*1.1);
        CriterionSimple criteria = new CriterionSimple(this, space, potential.getRange(), potential.getRange()*1.1);
        potentialMaster.setCriterion(potential, new CriterionTypesCombination(criteria, new AtomType[]{movable.getLeafType(), potentialSpecies.getLeafType()}));

        this.potentialMasterD.addPotential(potential, new AtomType[]{movable.getLeafType(), potentialSpecies.getLeafType()});
        potentialMasterD.setSpecies(new ISpecies []{potentialSpecies, movable});
        potentialMasterD.setRange(potential.getRange()*1.1);
        CriterionSimple criteria2 = new CriterionSimple(this, space, potential.getRange(), potential.getRange()*1.1);
        potentialMasterD.setCriterion(potential, new CriterionTypesCombination(criteria2, new AtomType[]{movable.getLeafType(), potentialSpecies.getLeafType()}));
        
    //ADATOM CREATION AND PLACEMENT
        // Sn
        
        IMolecule iMolecule = movable.makeMolecule();
        box.addMolecule(iMolecule);
        adAtomPos = iMolecule.getChildList().getAtom(0).getPosition();
        //adAtomPos = getSpace().makeVector();
        adAtomPos.setX(0, 9.8);
        adAtomPos.setX(1, 0.2);
        adAtomPos.setX(2, -1.0);
        Vector newBoxLength = space.makeVector();
        newBoxLength.E(box.getBoundary().getBoxSize());
        newBoxLength.setX(0, 2.0*adAtomPos.getX(0)+2.0);
        box.getBoundary().setBoxSize(newBoxLength);
        
        /**
        //Ag
        IAtom iAtom = agAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(agAdatom).addChildAtom(iAtom);
        iAtom).getPosition().setX(0, 7.5);
        iAtom).getPosition().setX(1, 0.9477016722828758);
        iAtom).getPosition().setX(2, 1.0709520701043456);
        */
        /**
        //Cu
        IAtom iAtom = cuAdatom.getMoleculeFactory().makeAtom();
        box.getAgent(cuAdatom).addChildAtom(iAtom);
        iAtom).getPosition().setX(0, 6.0);
        iAtom).getPosition().setX(1, 0.9477016722828758);
        iAtom).getPosition().setX(2, 1.0709520701043456);
        */
    }

    public static void main(String[] args) {

        final SimDimerMEAMadatom sim = new SimDimerMEAMadatom();
        Vector vect = sim.getSpace().makeVector();
        vect.setX(0, 9.8);
        vect.setX(1, -0.2);
        vect.setX(2, -0.2);

        sim.setMovableAtoms(100.0, vect);

        sim.setPotentialListAtoms();
        //sim.removeAtoms(2.9, vect);

        //sim.enableMolecularDynamics(5000);

        //sim.enableDimerSearch("0-MEAM", 2000, false, false);
        //sim.integratorDimer.setRotNum(0);
        sim.initializeConfiguration("0-MEAM_A_minimum");
        CalcVibrationalModes vib = new CalcVibrationalModes();
        vib.setup(sim.box, sim.potentialMasterD, sim.box.getMoleculeList(sim.movable), sim.getSpace());
        vib.actionPerformed();
        System.out.println(vib.getProductOfFrequencies());
        sim.initializeConfiguration("0-MEAM_saddle");
        vib.actionPerformed();
        System.out.println(vib.getProductOfFrequencies());
        //sim.enableMinimumSearch("0-MEAM", false);
        //sim.integratorDimerMin.initializeDimer();

        /*
        WriteConfiguration poswriter = new WriteConfiguration(sim.space);
        poswriter.setConfName("sns101-initial3");
        poswriter.setBox(sim.box);
        sim.integratorMD.addIntervalAction(poswriter);
        sim.integratorMD.setActionInterval(poswriter, 1);
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
        IntegratorListenerAction pumpListener = new IntegratorListenerAction(energyPump);
        pumpListener.setInterval(1);
        sim.integratorDimerMin.getEventManager().addListener(pumpListener);

        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1, sim.space, sim.getController());
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));
        simGraphic.add(plotPE);
        //sim.integratorMD.addIntervalAction(simGraphic.getPaintAction(sim.box));
        sim.integratorDimerMin.getEventManager().addListener(new IntegratorListenerAction(simGraphic.getPaintAction(sim.box)));

        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());

        colorScheme.setColor(sim.fixed.getLeafType(), java.awt.Color.gray);
        colorScheme.setColor(sim.movable.getLeafType(), java.awt.Color.red);
        colorScheme.setColor(sim.potentialSpecies.getLeafType(), java.awt.Color.PINK);

        simGraphic.makeAndDisplayFrame(APP_NAME);
    }
    
    public void setMovableAtoms(double distance, Vector center){
        //distance = distance*distance;
        Vector rij = space.makeVector();
        MoleculeArrayList movableList = new MoleculeArrayList();
        IMoleculeList loopSet = box.getMoleculeList();
        for (int i=0; i<loopSet.getMoleculeCount(); i++){
            rij.Ev1Mv2(center,loopSet.getMolecule(i).getChildList().getAtom(0).getPosition());
            if(rij.getX(0) > (box.getBoundary().getBoxSize().getX(0) - 3.0)){continue;}
            //box.getBoundary().nearestImage(rij);
            if(rij.squared() < distance){
               movableList.add(loopSet.getMolecule(i));
            }
        }
        for (int i=0; i<movableList.getMoleculeCount(); i++){
            IMolecule newMolecule = movable.makeMolecule();
            box.addMolecule(newMolecule);
            newMolecule.getChildList().getAtom(0).getPosition().E(movableList.getMolecule(i).getChildList().getAtom(0).getPosition());
            box.removeMolecule(movableList.getMolecule(i));
        }
        movableSet = box.getMoleculeList(movable);
    }
    
    public void setPotentialListAtoms(){
        MoleculeArrayList neighborList = new MoleculeArrayList();
        MoleculeArrayList fixedList = new MoleculeArrayList();
        IMoleculeList loopSet = box.getMoleculeList();
        for(int i=0; i<loopSet.getMoleculeCount(); i++){
            if(loopSet.getMolecule(i).getType()==movable){
                continue;
            }
            if(loopSet.getMolecule(i).getChildList().getAtom(0).getPosition().getX(0) < -4.0){continue;}
            boolean fixedFlag = true;
            for(int j=0; j<movableSet.getMoleculeCount(); j++){
                Vector dist = space.makeVector();
                dist.Ev1Mv2(loopSet.getMolecule(i).getChildList().getAtom(0).getPosition(),movableSet.getMolecule(j).getChildList().getAtom(0).getPosition());
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
            IMolecule newMolecule = potentialSpecies.makeMolecule();
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
        Vector rij = space.makeVector();

        IMoleculeList loopSet = box.getMoleculeList(movable);
        for (int i=0; i<loopSet.getMoleculeCount(); i++){
            rij.Ev1Mv2(center,loopSet.getMolecule(i).getChildList().getAtom(0).getPosition());
            box.getBoundary().nearestImage(rij);
            if(rij.squared() < distance){
               box.removeMolecule(loopSet.getMolecule(i));
            }
        }
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

    public void enableMolecularDynamics(long maxSteps){
        integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integratorMD.setTimeStep(0.001);
        integratorMD.setTemperature(Kelvin.UNIT.toSim(100));
        integratorMD.setThermostatInterval(100);
        integratorMD.setIsothermal(true);
        integratorMD.setBox(box);
        integratorMD.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        activityIntegrateMD = new ActivityIntegrate(integratorMD);
        getController().addAction(activityIntegrateMD);
        activityIntegrateMD.setMaxSteps(maxSteps);
    }

    public void enableDimerSearch(String fileName, long maxSteps, Boolean orthoSearch, Boolean fine){

        integratorDimer = new IntegratorDimerRT(this, potentialMasterD, new ISpecies[]{movable}, space);
        integratorDimer.setBox(box);
        integratorDimer.setOrtho(orthoSearch, false);
        if(fine){
            ConfigurationFile configFile = new ConfigurationFile(fileName+"_saddle");
            configFile.initializeCoordinates(box);

            integratorDimer.setFileName(fileName+"_fine");
            integratorDimer.deltaR = 0.0005;
            integratorDimer.dXl = 10E-5;
            integratorDimer.deltaXmax = 0.005;
            integratorDimer.dFsq = 0.0001*0.0001;
            integratorDimer.dFrot = 0.01;
        }
        integratorDimer.setFileName(fileName);
        integratorDimer.getEventManager().addListener(potentialMasterD.getNeighborManager(box));
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
        integratorDimer.setActivityIntegrate(activityIntegrateDimer);
        getController().addAction(activityIntegrateDimer);
        activityIntegrateDimer.setMaxSteps(maxSteps);
    }

    public void enableMinimumSearch(String fileName, Boolean normalDir){

        integratorDimerMin = new IntegratorDimerMin(this, potentialMasterD, new ISpecies[]{movable}, normalDir, space);
        integratorDimerMin.setBox(box);
        integratorDimerMin.setFileName(fileName);
        integratorDimerMin.getEventManager().addListener(potentialMasterD.getNeighborManager(box));
        activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
        integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
        getController().addAction(activityIntegrateMin);
    }

}
