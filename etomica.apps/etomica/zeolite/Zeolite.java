package etomica.zeolite;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationSequential;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionSimple;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.P1HardBoundary;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;
import etomica.config.ConfigurationFile;
import etomica.potential.P2LennardJones;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.BoundaryRectangularPeriodic;


public class Zeolite {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		Default defaults = new Default();
        defaults.doSleep = false;
        defaults.ignoreOverlap = true;
        ZeoliteSimulation[] sim = new ZeoliteSimulation[1];
        for(int i=0;i<sim.length;i++){
        	sim[i]=new ZeoliteSimulation(defaults);
        	sim[i].getController().actionPerformed();
        	System.out.println("Simulations Complete");
        	String file = sim[i].getFileName();
        	Converter(file,sim[i].getMethane(),i+11);
        }
        
        
        
        System.out.println("end");
	}
	public static void Converter(String inputFile,int meth,int index) {
		// TODO Auto-generated method stub
		String outputFile = inputFile+"_"+index+"_Result.txt";
		MSDProcessor proc = new MSDProcessor(Space3D.getInstance(),inputFile,outputFile);
		
		proc.setDeltaTmax(1);
		proc.setMethane(meth);
		proc.fillArrays();
		System.out.println("Converter done");
	}

}
