// Source file generated by Etomica

package etomica.simulation.prototypes;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.box.Box;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorHard;
import etomica.lattice.LatticeCubicFcc;
import etomica.potential.P2HardSphere;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;

public class HSMD3DNoNbr extends Simulation {

    private static final long serialVersionUID = 1L;
    public IBox box;
    public IntegratorHard integrator;
    public SpeciesSpheresMono species;
    public P2HardSphere potential;
    
    public HSMD3DNoNbr() {
        super(Space3D.getInstance(), true);
        IPotentialMaster potentialMaster = new PotentialMaster(space);

        int numAtoms = 256;
        double sigma = 1.0;
        double l = 14.4573*Math.pow((numAtoms/2020.0),1.0/3.0);

        integrator = new IntegratorHard(this, potentialMaster, space);
        integrator.setIsothermal(false);
        integrator.setTimeStep(0.01);

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(integrator);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this, space);
        getSpeciesManager().addSpecies(species);
        potential = new P2HardSphere(space, sigma, false);
        potentialMaster.addPotential(potential,new IAtomTypeLeaf[]{species.getLeafType(),species.getLeafType()});

        box = new Box(this, space);
        addBox(box);
        box.setNMolecules(species, numAtoms);
        box.setDimensions(space.makeVector(new double[]{l,l,l}));
//        box.setBoundary(new BoundaryTruncatedOctahedron(space));
        integrator.setBox(box);
        integrator.addIntervalAction(new BoxImposePbc(box, space));
        new ConfigurationLattice(new LatticeCubicFcc(), space).initializeCoordinates(box);
        
        //ColorSchemeByType.setColor(speciesSpheres0, java.awt.Color.blue);

 //       MeterPressureHard meterPressure = new MeterPressureHard(integrator);
 //       DataAccumulator accumulatorManager = new DataAccumulator(meterPressure);
        // 	DisplayTextBox box = new DisplayBox();
        // 	box.setDatumSource(meterPressure);
 //       box.setDensity(0.7);
    } //end of constructor

    public static void main( String[] args )
    {
    	String filename = "test.bin";
		
    	try
    	{
    	    FileOutputStream fos = null;
    	    ObjectOutputStream out = null;
    	    HSMD3DNoNbr simulation = new HSMD3DNoNbr();
    	    fos = new FileOutputStream( filename);
			out = new ObjectOutputStream(fos);
			out.writeObject( simulation );
			out.close();
			fos.close();
			System.out.println( "Serialization of class HSMD3DNoNbr succeeded.");
    	}
    	catch(IOException ex)
    	{
    	    System.err.println( "Exception:" + ex.getMessage() );
    	    ex.printStackTrace();
    	}
    	
    	// Serialize back
    	ISimulation simulation = null;
    	try
    	{
    	    FileInputStream fis = null;
    	    fis = new FileInputStream(filename);
    	    ObjectInputStream in = new ObjectInputStream(fis);
    	    simulation = (ISimulation) in.readObject();
    	    in.close();
    	    fis.close();
    	    
    	    System.out.println( "DeSerialization of class HSMD3DNoNbr succeeded.");

    	}
    	catch( Exception ex ) {
    	    System.err.println( "Could not read simulation from file " + filename + ". Cause: " + ex.getMessage() );
    	    ex.printStackTrace();
    	}
		
	    // go daddy
	    simulation.getController().actionPerformed();
	    System.out.println( "Simulation run ok");
		
    }
}//end of class
