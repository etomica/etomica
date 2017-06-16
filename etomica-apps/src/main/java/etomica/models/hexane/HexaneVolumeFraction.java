/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.hexane;

/**
 * Class to calculate the volume of a single hexane molecule
 */
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.atom.Atom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.space3d.Space3D;

public class HexaneVolumeFraction extends Simulation {

    public ActivityIntegrate activityIntegrate;
    public IntegratorMC integrator;

    public Box box;
    public BoundaryRectangularPeriodic bdry;
    public MCMoveMolecule moveMolecule;
       
    public HexaneVolumeFraction(Space _space) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(_space);
        PotentialMaster potentialMaster = new PotentialMaster();
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = 6;


        SpeciesHexane species = new SpeciesHexane(_space);
        addSpecies(species);

        bdry = new BoundaryRectangularPeriodic(_space);
        box = new Box(bdry, _space);
        addBox(box);        
        box.getBoundary().setBoxSize(space.makeVector(new double[] {3.8, 3.8, 3.8}));
        box.setNMolecules(species, 1);

        integrator = new IntegratorMC(potentialMaster, getRandom(), 1.0);

        activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setMaxSteps(2000000);
        getController().addAction(activityIntegrate);
  
        integrator.setBox(box);
        

       
    }
    /**
     * @param args
     */
    public static void main(String[] args) {
        boolean graphic = false;
        double numberOfTests = 1000000000;
        
        HexaneVolumeFraction sim = new HexaneVolumeFraction(Space3D.getInstance());
        RandomPositionSource positionSource = new RandomPositionSourceRectangular(sim.space, sim.getRandom());
        positionSource.setBox(sim.box);
        
        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space, sim.getController());
            simGraphic.makeAndDisplayFrame();
        } else {
            long overlaps = 0;
            long time = System.currentTimeMillis();
            System.out.println(time);
            long time1;
            long time2;

            Vector temp = sim.space.makeVector();
            Vector rand = sim.space.makeVector();
            AtomIteratorLeafAtoms ail = new AtomIteratorLeafAtoms(sim.box);
            ail.reset();
            IAtom atom = new Atom(sim.getSpace());
            
            time1 = System.currentTimeMillis();
            for(int count = 0; count < numberOfTests; count++){
                rand = positionSource.randomPosition();
                ail.reset();

//                //FIRST METHOD OF LOOPING
//                while((atom = (AtomLeaf)ail.nextAtom()) != null){
//                    temp.E(((AtomLeaf)atom).getPosition());
//                    temp.ME(rand);
//                    double length = Math.sqrt(temp.dot(temp));
//                    if(length <= 0.5){
//                        overlaps += 1;
//                        break;
//                    }
//                }
//                
//                //SECOND METHOD OF LOOPING
//                for(atom = (AtomLeaf)ail.nextAtom(); atom != null; atom = (AtomLeaf)ail.nextAtom()){
//                    temp.E(((AtomLeaf)atom).getPosition());
//                    temp.ME(rand);
//                    double length = Math.sqrt(temp.dot(temp));
//                    if(length <= 0.5){
//                        overlaps++;
//                        break;
//                    }
//                }
                                
                //THIRD METHOD OF LOOPING
                IAtomList list = sim.box.getLeafList();
                for(int i = 0; i < list.getAtomCount(); i++){
                    atom = list.getAtom(i);
                    temp.E(atom.getPosition());
                    temp.ME(rand);
                    double length = Math.sqrt(temp.dot(temp));
                    if(length <= 0.5){
                        overlaps++;
                        break;
                    }
                }
                
            }

            time2 = System.currentTimeMillis();
            System.out.println("start:  " + time);
            System.out.println("simulation:  " + time1);
            System.out.println("data collection:  " + time2);
            System.out.println("overlaps:  " + overlaps);
            System.out.println("total tries:  " + numberOfTests);
            double td = (double)overlaps/ numberOfTests;
            System.out.println("percentage in molecule:  " + td);
            double vol = sim.bdry.volume();
            System.out.println("volume of box:  " + vol);
            double another = td * vol;
            System.out.println("volume of molecule:  " + another);
        
        }
    }
}
