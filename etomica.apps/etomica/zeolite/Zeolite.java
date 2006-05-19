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
		// Converts the simulation data into something I can use
		String inputFile = "32_1000000_0.00611_2000_WCA";
		//Cut data file down to multiple parts
		int cut = 5;
		DataCutter data = new DataCutter(inputFile,cut);
		
		for(int i=0;i<cut;i++){
			String file = inputFile+"_"+i;
			Converter(file,32);
			System.out.println("File "+i+" converted");
		}
		
		//Converter(inputFile,32,2);     
        System.out.println("Finished");
	}
	public static void Converter(String inputFile,int meth) {
		// TODO Auto-generated method stub
		String outputFile = inputFile+"_Result.txt";
		MSDProcessor proc = new MSDProcessor(Space3D.getInstance(),inputFile,outputFile);
		
		//proc.setDeltaTmax(1);
		proc.setMethane(meth);
		proc.fillArrays();
		System.out.println("Converter done");
	}

}
