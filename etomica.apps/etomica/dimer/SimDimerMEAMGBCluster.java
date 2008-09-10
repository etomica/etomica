package etomica.dimer;

import etomica.action.WriteConfiguration;
import etomica.action.XYZWriter;
import etomica.api.IVector;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.DisplayPlot;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.util.HistoryCollapsingAverage;

/**
 * Simulation using Henkelman's Dimer method to find a saddle point for
 * an adatom of Sn on a surface, modeled with MEAM.
 * 
 * @author msellers
 *
 */

public class SimDimerMEAMGBCluster extends Simulation{
	

	public SimDimerMEAMGBCluster(Space space) {
		super(space);
	}

	public static void main(String[] args){
		
		String fileName = args[0];
        //int mdSteps = 10;//Integer.parseInt(args[1]);
        //int h = Integer.parseInt(args[1]);
        //int k = Integer.parseInt(args[2]);
        //int l = Integer.parseInt(args[3]);
        
        //int x = Integer.parseInt(args[4]);
        //int y = Integer.parseInt(args[5]);
        //int z = Integer.parseInt(args[6]);
        
        int num = Integer.parseInt(args[1]);
        
        final String APP_NAME = "SimDimerMEAMGBCluster";
        
        final SimDimerMEAMGB sim = new SimDimerMEAMGB(new int[] {2,1,0}, new int[] {2,6,12});
        
        sim.initializeConfiguration("sngb210-2612-md");
        
        
        IVector dimerCenter = sim.getSpace().makeVector();
        dimerCenter.setX(0, sim.box.getBoundary().getDimensions().x(0)/2.0);
        dimerCenter.setX(1, 1.0);
        dimerCenter.setX(2, 0.0);
         
        if(sim.millerPlane[2] == 0){
            dimerCenter.setX(1, sim.box.getBoundary().getDimensions().x(1)/2.0);
            dimerCenter.setX(0, 1.0);
            dimerCenter.setX(2, 0.0);
        }
        
        sim.setMovableAtoms(6.0, dimerCenter);
        sim.setMovableAtomsList();
        
        /*
        sim.initializeConfiguration(fileName+"_saddle");
        sim.calculateVibrationalModes(fileName+"_saddle");
        
        sim.initializeConfiguration(fileName+"_A_minimum");
        sim.calculateVibrationalModes(fileName+"_A_minimum");
        
        sim.initializeConfiguration(fileName+"_B_minimum");
        sim.calculateVibrationalModes(fileName+"_B_minimum");
        */      
        sim.initializeConfiguration(fileName+"_saddle");
        
        //sim.enableMolecularDynamics(10000);
        
        //sim.enableDimerSearch(fileName, 2500, false, false);
        //sim.integratorDimer.setRotNum(num);
        
        sim.enableMinimumSearch(fileName, true);
        
        XYZWriter xyzwriter = new XYZWriter(sim.box);
        xyzwriter.setFileName(fileName+"-B_minumum.xyz");
        xyzwriter.setIsAppend(true);
        sim.integratorDimerMin.addIntervalAction(xyzwriter);
        sim.integratorDimerMin.setActionInterval(xyzwriter, 5);
        
        /*
        WriteConfiguration writer = new WriteConfiguration(sim.getSpace());
        writer.setBox(sim.box);
        writer.setConfName(fileName+"-md");
        sim.integratorMD.addIntervalAction(writer);
        sim.integratorMD.setActionInterval(writer, 10000);
		*/
        
        //MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
        //energyMeter.setBox(sim.box);
         
        sim.getController().actionPerformed();

    }

}
