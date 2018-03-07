/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.action.BoxImposePbc;
import etomica.action.BoxInflate;
import etomica.action.CalcVibrationalModes;
import etomica.action.WriteConfiguration;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationFile;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.integrator.IntegratorListenerAction;
import etomica.math.numerical.CalcGradientDifferentiable;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeArrayList;
import etomica.potential.P2LennardJones;
import etomica.potential.PotentialMaster;
import etomica.potential.PotentialMasterMonatomic;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom on a surface, modeled with LJ.
 * 
 * @author msellers
 *
 */

public class SimDimerLJadatom extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "DimerLJadatom";
    public final PotentialMaster potentialMaster;
    public IntegratorVelocityVerlet integratorMD;
    public IntegratorDimerRT integratorDimer;
    public IntegratorDimerMin integratorDimerMin;
    public Box box;
    public Vector[] saddle, normal;
    public SpeciesSpheresMono fixed, movable;
//    public P2LennardJones potential;
    public ActivityIntegrate activityIntegrateMD, activityIntegrateDimer, activityIntegrateMin;
    public CalcGradientDifferentiable calcGradientDifferentiable;
    public CalcVibrationalModes calcVibrationalModes;
    public double [][] dForces;
    public int [] d, modeSigns;
    public double [] positions;
    public double [] lambdas, frequencies;
    public IMoleculeList movableSet;
    public Vector adAtomPos;
    public Boolean saddleFine, calcModes, minSearch, normalDir;
    

    public SimDimerLJadatom() {
        super(Space3D.getInstance());
        potentialMaster = new PotentialMasterMonatomic(this);

        //SIMULATION BOX
        box = this.makeBox(new BoundaryRectangularSlit(0, 5, space));

        //SPECIES
        double sigma = 1.0;
        fixed = new SpeciesSpheresMono(space, new ElementSimple("A", Double.POSITIVE_INFINITY));
        fixed.setIsDynamic(true);
        movable = new SpeciesSpheresMono(this, space);
        movable.setIsDynamic(true);
        addSpecies(fixed);
        addSpecies(movable);

        // Must be in same order as the respective species is added to SpeciesManager
        box.setNMolecules(fixed, 256);

        BoxInflate inflater = new BoxInflate(box, space);
        inflater.setTargetDensity(1);
        inflater.actionPerformed();

//    	potential = new P2LennardJones(space, sigma, 1.0);
//		potentialMaster.addPotential(, new IAtomTypeLeaf[]{fixed.getLeafType(), fixed.getLeafType()});
        potentialMaster.addPotential(new P2LennardJones(space, sigma, 1.0), new AtomType[]{movable.getLeafType(), fixed.getLeafType()});
        potentialMaster.addPotential(new P2LennardJones(space, sigma, 1.0), new AtomType[]{movable.getLeafType(), movable.getLeafType()});

        //CRYSTAL
        Configuration config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);

        //ADATOM CREATION AND PLACEMENT

        IMolecule iMolecule = movable.makeMolecule();
        box.addMolecule(iMolecule);
        adAtomPos = iMolecule.getChildList().get(0).getPosition();
        //adAtomPos = getSpace().makeVector();
        adAtomPos.setX(0, 3.5);
        adAtomPos.setX(1, -0.30);
        adAtomPos.setX(2, -0.30);
        Vector newBoxLength = space.makeVector();
        newBoxLength.E(box.getBoundary().getBoxSize());
        newBoxLength.setX(0, 2.0 * adAtomPos.getX(0) + 2.0);
        box.getBoundary().setBoxSize(newBoxLength);

    }

    public static void main(String[] args) {

        final SimDimerLJadatom sim = new SimDimerLJadatom();
        Vector vect = sim.getSpace().makeVector();
        vect.setX(0, 3.5);
        vect.setX(1, 0.0);
        vect.setX(2, 0.0);

        sim.setMovableAtoms(2.0, vect);

        //sim.initializeConfiguration("0");

        //sim.enableDimerSearch("0", 800, false, false);
        //sim.randomizePositions();
        //sim.integratorDimer.setRotNum(0);
        sim.enableMinimumSearch("s_1", false);


        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME, 1);
        simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

        sim.integratorDimerMin.getEventManager().addListener(new IntegratorListenerAction(simGraphic.getPaintAction(sim.box)));

        ColorSchemeByType colorScheme = ((ColorSchemeByType) ((DisplayBox) simGraphic.displayList().getFirst()).getColorScheme());

        colorScheme.setColor(sim.fixed.getLeafType(), java.awt.Color.gray);
        colorScheme.setColor(sim.movable.getLeafType(), java.awt.Color.red);

        simGraphic.makeAndDisplayFrame(APP_NAME);
    }

    public void setMovableAtoms(double distance, Vector center){
        //distance = distance*distance;
        Vector rij = space.makeVector();
        MoleculeArrayList movableList = new MoleculeArrayList();
        IMoleculeList loopSet = box.getMoleculeList();
        for (int i = 0; i<loopSet.size(); i++){
            rij.Ev1Mv2(center,loopSet.get(i).getChildList().get(0).getPosition());
            if(rij.getX(0) > (box.getBoundary().getBoxSize().getX(0) - 3.0)){continue;}
            //box.getBoundary().nearestImage(rij);
            if(rij.getX(0)< distance){
               movableList.add(loopSet.get(i));
            }
        }
        for (int i = 0; i<movableList.size(); i++){
            IMolecule newMolecule = movable.makeMolecule();
            box.addMolecule(newMolecule);
            newMolecule.getChildList().get(0).getPosition().E(movableList.get(i).getChildList().get(0).getPosition());
            box.removeMolecule(movableList.get(i));
        }
        movableSet = box.getMoleculeList(movable);
    }
    
    //Must be run after setMovableAtoms
    public void removeAtoms(double distance, Vector center){
        distance = distance*distance;
        Vector rij = space.makeVector();

        IMoleculeList loopSet = box.getMoleculeList(movable);
        for (int i = 0; i<loopSet.size(); i++){
            rij.Ev1Mv2(center,loopSet.get(i).getChildList().get(0).getPosition());
            box.getBoundary().nearestImage(rij);
            if(rij.squared() < distance){
               box.removeMolecule(loopSet.get(i));
            }
        }
    }
    
    public void initializeConfiguration(String fileName){
        ConfigurationFile config = new ConfigurationFile(fileName);
        config.initializeCoordinates(box);
    }

    public void generateConfigs(String fileName, double percentd) {

        Vector workVector = space.makeVector();
        Vector[] currentPos = new Vector[movableSet.size()];
        for(int i=0; i<currentPos.length; i++){
            currentPos[i] = space.makeVector();
            currentPos[i].E(movableSet.get(i).getChildList().get(0).getPosition());
        }

        //Create multiple configurations
        for(int m=0; m<50; m++){
            WriteConfiguration genConfig = new WriteConfiguration(space);
            genConfig.setBox(box);
            genConfig.setConfName(fileName+"_config_"+m);
            //Displaces atom's by at most +/-0.03 in each coordinate
            for(int i = 0; i<movableSet.size(); i++){
                Vector atomPosition = movableSet.get(i).getChildList().get(0).getPosition();
                for(int j=0; j<3; j++){
                    workVector.setX(j,percentd*random.nextGaussian());
                }
                atomPosition.Ev1Pv2(currentPos[i],workVector);
            }
            genConfig.actionPerformed();
        }
    }

    public void enableMolecularDynamics(long maxSteps) {
        integratorMD = new IntegratorVelocityVerlet(this, potentialMaster, box);
        integratorMD.setTimeStep(0.01);
        integratorMD.setTemperature(0.1);
        integratorMD.setThermostatInterval(100);
        integratorMD.setIsothermal(true);
        activityIntegrateMD = new ActivityIntegrate(integratorMD);
        BoxImposePbc imposePbc = new BoxImposePbc(box, space);
        integratorMD.getEventManager().addListener(new IntegratorListenerAction(imposePbc));
        getController().addAction(activityIntegrateMD);
        activityIntegrateMD.setMaxSteps(maxSteps);
    }

    public void enableDimerSearch(String fileName, long maxSteps, Boolean orthoSearch, Boolean fine) {

        integratorDimer = new IntegratorDimerRT(this, potentialMaster, new ISpecies[]{movable}, box);
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
        activityIntegrateDimer = new ActivityIntegrate(integratorDimer);
        integratorDimer.setActivityIntegrate(activityIntegrateDimer);
        getController().addAction(activityIntegrateDimer);
        activityIntegrateDimer.setMaxSteps(maxSteps);
    }

    public void enableMinimumSearch(String fileName, Boolean normalDir) {

        integratorDimerMin = new IntegratorDimerMin(this, potentialMaster, new ISpecies[]{movable}, normalDir, box);
        integratorDimerMin.setFileName(fileName);
        activityIntegrateMin = new ActivityIntegrate(integratorDimerMin);
        integratorDimerMin.setActivityIntegrate(activityIntegrateMin);
        getController().addAction(activityIntegrateMin);
    }
    
    public void randomizePositions(){
        Vector workVector = space.makeVector();
        IMoleculeList loopSet3 = box.getMoleculeList(movable);
        Vector[] currentPos = new Vector[loopSet3.size()];
        double offset = 0;
        for(int i=0; i<currentPos.length; i++){
            currentPos[i] = space.makeVector();
            currentPos[i] = (loopSet3.get(i).getChildList().get(0).getPosition());
            for(int j=0; j<3; j++){
                offset = random.nextGaussian()/10.0;
                if(Math.abs(offset)>0.1){offset=0.1;}
                workVector.setX(j,offset);
            }
            currentPos[i].PE(workVector);
        }
    }

}
