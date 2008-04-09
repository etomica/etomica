package etomica.models.hexane;

/**
 * Class to calculate the volume of a single hexane molecule
 */
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
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

    public IBox box;
    public BoundaryRectangularPeriodic bdry;
    public MCMoveMolecule moveMolecule;
       
    public HexaneVolumeFraction(Space _space) {
        //super(space, false, new PotentialMasterNbr(space, 12.0));
//        super(space, true, new PotentialMasterList(space, 12.0));
        super(_space, false);
        PotentialMaster potentialMaster = new PotentialMaster(_space);
        int chainLength = 6;
        //One molecule per cell
        int numAtoms = 6;


        SpeciesHexane species = new SpeciesHexane(this, _space);
        getSpeciesManager().addSpecies(species);

        bdry = new BoundaryRectangularPeriodic(getRandom(), _space);
        box = new Box(bdry, _space);
        addBox(box);        
        box.setDimensions(space.makeVector(new double[] {3.8, 3.8, 3.8}));
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
        
        if (graphic) {
            SimulationGraphic simGraphic = new SimulationGraphic(sim, sim.space);
            simGraphic.makeAndDisplayFrame();
        } else {
            long overlaps = 0;
            long time = System.currentTimeMillis();
            System.out.println(time);
            long time1;
            long time2;

            IVector temp = sim.space.makeVector();
            IVector rand = sim.space.makeVector();
            AtomIteratorLeafAtoms ail = new AtomIteratorLeafAtoms(sim.box);
            ail.reset();
            AtomLeaf atom = new AtomLeaf(sim.getSpace());
            
            time1 = System.currentTimeMillis();
            for(int count = 0; count < numberOfTests; count++){
                rand = ((BoundaryRectangularPeriodic)sim.box.getBoundary()).randomPosition();
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
                AtomArrayList list = (AtomArrayList)sim.box.getLeafList();
                for(int i = 0; i < list.getAtomCount(); i++){
                    atom = (AtomLeaf)list.getAtom(i);
                    temp.E(((AtomLeaf)atom).getPosition());
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
            double td = (double)overlaps/(double)numberOfTests;
            System.out.println("percentage in molecule:  " + td);
            double vol = ((BoundaryRectangularPeriodic)sim.bdry).volume();
            System.out.println("volume of box:  " + vol);
            double another = td * vol;
            System.out.println("volume of molecule:  " + another);
        
        }
    }
}
