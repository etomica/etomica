/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.dimer;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorListenerAction;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.meam.ParameterSetMEAM;
import etomica.meam.PotentialMEAM;
import etomica.molecule.IMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;

/**
 * Simulation class using IntegratorEnergyMap.  Rasterizes a surface with an atom
 * to achieve values of its energy.
 * 
 * @author msellers
 *
 */

public class EnergyMap extends Simulation{

    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM Md3D";
    public final PotentialMaster potentialMaster;
    public IntegratorEnergyMap integratorMAP;
    public Box box;
    public SpeciesSpheresMono sn, snFix, snAdatom, cu, cuFix, cuAdatom, movable;
    public PotentialMEAM potential;
    public ActivityIntegrate activityIntegrateMAP;
    

    
    public EnergyMap(double height, String fileTail) {
        super(Space3D.getInstance());

        // Sn
        Tin tinFixed = new Tin("SnFix", Double.POSITIVE_INFINITY);

        snFix = new SpeciesSpheresMono(space, tinFixed);
        snFix.setIsDynamic(true);
        sn = new SpeciesSpheresMono(space, Tin.INSTANCE);
        sn.setIsDynamic(true);
        snAdatom = new SpeciesSpheresMono(space, Tin.INSTANCE);
        snAdatom.setIsDynamic(true);
        movable = new SpeciesSpheresMono(space, Tin.INSTANCE);
        movable.setIsDynamic(true);

        addSpecies(snFix);
        addSpecies(sn);
        addSpecies(snAdatom);
        addSpecies(movable);

        /**
         // Cu
         Copper copperFixed = new Copper("CuFix", Double.POSITIVE_INFINITY);

         cuFix = new SpeciesSpheresMono(this, copperFixed);
         cu = new SpeciesSpheresMono(this, Copper.INSTANCE);
         cuAdatom = new SpeciesSpheresMono(this, Copper.INSTANCE);
         movable = new SpeciesSpheresMono(this, Copper.INSTANCE);

         getSpeciesManager().addSpecies(cuFix);
         getSpeciesManager().addSpecies(cu);
         getSpeciesManager().addSpecies(cuAdatom);
         getSpeciesManager().addSpecies(movable);

         */

        potentialMaster = new PotentialMaster();

        box = this.makeBox(new BoundaryRectangularSlit(0, 5, space));

        // Sn
        box.setNMolecules(snFix, 72);
        box.setNMolecules(sn, 144);
        box.setNMolecules(snAdatom, 0);

        potential = new PotentialMEAM(space);

        potential.setParameters(snFix.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(sn.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(snAdatom.getLeafType(), ParameterSetMEAM.Sn);
        potential.setParameters(movable.getLeafType(), ParameterSetMEAM.Sn);

        this.potentialMaster.addPotential(potential, new AtomType[]{sn.getLeafType(), snFix.getLeafType(), snAdatom.getLeafType(), movable.getLeafType()});

        /**
         // Cu
         box.setNMolecules(cuFix, 64);
         box.setNMolecules(cu, 128);
         box.setNMolecules(cuAdatom, 0);

         potential = new PotentialMEAM(space);

         potential.setParameters(cuFix, ParameterSetMEAM.Cu);
         potential.setParameters(cu, ParameterSetMEAM.Cu);
         potential.setParameters(cuAdatom, ParameterSetMEAM.Cu);
         potential.setParameters(movable, ParameterSetMEAM.Cu);

         this.potentialMaster.addPotential(potential, new Species[]{cu, cuFix, cuAdatom, movable});
         */

        // Sn box

        box.getBoundary().setBoxSize(new Vector3D(5.8314 * 3, 5.8314 * 3, 3.1815 * 6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisBetaSnA5());

        /**
         // Cu box
         box.setDimensions(new Vector3D(3.6148*3, 3.6148*4, 3.6148*4));
         PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
         BravaisLatticeCrystal crystal = new BravaisLatticeCrystal(primitive, new BasisCubicFcc());

         */


        Configuration config = new ConfigurationLattice(crystal, space);
        config.initializeCoordinates(box);

        // Sn
        IMolecule adMolecule = snAdatom.makeMolecule();
        box.addMolecule(adMolecule);
        IAtom adAtom = adMolecule.getChildList().get(0);

        adAtom.getPosition().setX(0, height);
        adAtom.getPosition().setX(1, -4.37);
        adAtom.getPosition().setX(2, -1.6);

        /**
         // Cu
         IAtom iAtom = cuAdatom.getMoleculeFactory().makeAtom();
         box.getAgent(cuAdatom).addChildAtom(iAtom);

         ((IAtomPositioned)iAtom).getPosition().setX(0, 6.3);
         ((IAtomPositioned)iAtom).getPosition().setX(1, -0.8);
         ((IAtomPositioned)iAtom).getPosition().setX(2, -1.2);
         */

        /**

         IVector rij = space.makeVector();
         AtomArrayList movableList = new AtomArrayList();
         AtomSet loopSet = box.getMoleculeList(sn);

         for (int i=0; i<loopSet.getAtomCount(); i++){
         rij.Ev1Mv2(((IAtomPositioned)iAtom).getPosition(),((IAtomPositioned)loopSet.getAtom(i)).getPosition());

         if(rij.squared()<32.49){
         movableList.add(loopSet.getAtom(i));
         }
         }

         for (int i=0; i<movableList.getAtomCount(); i++){
         ((IAtomPositioned)box.addNewMolecule(movable)).getPosition().E(((IAtomPositioned)movableList.getAtom(i)).getPosition());
         box.removeMolecule(movableList.getAtom(i));
         }

         */

        integratorMAP = new IntegratorEnergyMap(this, potentialMaster, adAtom, fileTail, box);
        activityIntegrateMAP = new ActivityIntegrate(integratorMAP);

        getController().addAction(activityIntegrateMAP);

    }

    public static void main(String[] args){
        double height1 = 10.0;
        String fileTail1 = ""+height1;
    	final String APP_NAME = "EnergyMapMEAMadatomSn";
    	final EnergyMap sim = new EnergyMap(height1, fileTail1);
    	
    	sim.activityIntegrateMAP.setMaxSteps(1); 
    	SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
    	simGraphic.getController().getReinitButton().setPostAction(simGraphic.getPaintAction(sim.box));

    	sim.integratorMAP.getEventManager().addListener(new IntegratorListenerAction(simGraphic.getPaintAction(sim.box)));
    	
    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());

    	//Sn
    
    	colorScheme.setColor(sim.sn.getLeafType(),java.awt.Color.gray);
        colorScheme.setColor(sim.snFix.getLeafType(),java.awt.Color.blue);
        colorScheme.setColor(sim.snAdatom.getLeafType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getLeafType(),java.awt.Color.PINK);
      
        /**
        
        //Cu
    
        colorScheme.setColor(sim.cu.getMoleculeType(),java.awt.Color.yellow);
        colorScheme.setColor(sim.cuFix.getMoleculeType(),java.awt.Color.cyan);
        colorScheme.setColor(sim.cuAdatom.getMoleculeType(),java.awt.Color.red);
        colorScheme.setColor(sim.movable.getMoleculeType(),java.awt.Color.PINK);
        
         */
    	
    	simGraphic.makeAndDisplayFrame(APP_NAME);
    	
    }

}
