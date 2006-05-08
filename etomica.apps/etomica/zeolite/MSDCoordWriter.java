package etomica.zeolite;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorMolecule;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.Species;
import etomica.units.*;
import etomica.data.meter.MeterTemperature;
/* =====SUMMARY======
 * At each 'writeInterval', which corresponds to a certain number of simulation steps,
 * every atom's absolute distance traveled is written to a file of the name 'fileName'.
 * Subsequent 'writeInverval' outputs append the file to create repeating blocks of the 
 * simulation's atoms with columns of X, Y, and Z.
 * ==================
 * 
 * The class is ordered as follows.  After the atoms have moved, but before periodic
 * boundary conditions are imposed, the main class's 'intervalAction' method gets atom
 * coordinates and stores them in 'atomOldCoord'.  The 'intervalCount', initially set to
 * 'writeInterval', is predecremented each simulation step.  Each time it does NOT equal
 * zero, the class continues and 'intervalAction' of the subclass 'AfterPBC' is called.  
 * The method fills 'atomPBIarray' using the pre-periodic boundary* 'atomOldCoord' and 
 * post-periodic boundary** 'workVector'. Each simulation step repeats this process.  When
 * 'intervalCount' IS zero, 'atomPBIarray' and the atoms' current positions are used as
 * stated below to create the file output.   
 * 
 *   
 * -----'atomPBIarray'----
 * Called in AfterPBC's 'intervalAction', this is an integer array containing a running
 * total of the number of box lengths traversed by the atoms during the simulation.  At 
 * a file write, 'atomPBIarray' is used as a multiplier of the primary box length.  
 * This quantity is added to the atoms' current coordinates to generate an atom's 
 * "absolute displacement" value.  
 * -----------------------
 * 
 *    * and ** : notice the getPriority() methods and their returns
 */


public class MSDCoordWriter implements IntegratorIntervalListener,
		IntegratorNonintervalListener {
	
	private int nAtomsMeth;
	private int totalAtoms;
	private Integrator integrate;
	private MeterTemperature meter;
	private double offset;
	
	public MSDCoordWriter(Space space, String fileName,Species[] species){
		// Creates an instance of subclass AfterPBC
		afterPBCinstance = new AfterPBC(space,species);
		iterator = new AtomIteratorMolecule(species);
		this.fileName = fileName;
		setWriteInterval(1);
		nAtomsMeth = 0;
		counter = 0;
	}
	
	public void setPhase(Phase phase){
		
		this.phase = phase;
		iterator.setPhase(phase);
		afterPBCinstance.setPhase(phase);
		totalAtoms = phase.atomCount();
		meter = new MeterTemperature();
		meter.setPhase(phase);
	}
	
	public void setNatoms(int natoms){
		nAtomsMeth = natoms;
		offset = totalAtoms/nAtomsMeth;
	}
	public void setIntegrator(Integrator integrator){
		integrate = integrator;
		integrator.addListener(this);
		integrator.addListener(afterPBCinstance);
	}
	
	// Methods involved with file creation/closing
	public void openFile(){
		try { 
			fileWriter = new FileWriter(fileName, false);
			fileWriter.write(iterator.size()+"\n");
		}
	    catch(IOException e) {
            System.err.println("Cannot open a file, caught IOException: " + e.getMessage());
        }
	}
	
	public void closeFile(){
        try {
            fileWriter.close();
            fileWriter = null;
        }
        catch(IOException e) {
            System.err.println("Cannot close a file, caught IOException: " + e.getMessage());
        }
    }
	
	/**
	 * Equates intervalAction method variable to variable passed in from simulation class
	 * @param writeInterval
	 */
	public void setWriteInterval(int writeInterval){
		this.writeInterval = writeInterval;
		intervalCount = writeInterval;
	}
	
	public void intervalAction(IntegratorIntervalEvent evt) {
		//System.out.println(counter);
		//counter++;
		afterPBCinstance.updateAtomOldCoord();
		if (--intervalCount == 0){
			counter++;
			System.out.println((counter*100*writeInterval/2000000.0)+"%");
			Vector phasedim = phase.getBoundary().getDimensions();
			// Gets atomPBIarray from AfterPBC subclass, through the subclass instance
			int [][] atomPBIarray = afterPBCinstance.getAtomPBIarray();
			double temp = Kelvin.UNIT.fromSim(meter.getDataAsScalar()*offset);
			try {
				fileWriter.write(temp+"\n");
				iterator.reset();
				int i=0;
				while (iterator.hasNext()){
					Vector atomPosition = ((AtomLeaf)iterator.nextAtom()).coord.position();
					for (int j=0;j < phasedim.D();j++){
						double actualDistance;
							
						// Total distance traveled between file writes is computed
						actualDistance = atomPBIarray[i][j] * phasedim.x(j) + atomPosition.x(j);
						fileWriter.write(""+actualDistance);
						if(j!=phasedim.D()-1){
							fileWriter.write("\t");
						}
						
					}
					fileWriter.write("\n");
					i++;
				}
			}
			catch (IOException e) {
	            throw new RuntimeException(e);
	        }
			// Variable is reset after a file write
			intervalCount = writeInterval;
		}
	}
	
	// *
	public int getPriority() {
		return 50;
	}

	public void nonintervalAction(IntegratorNonintervalEvent evt) {

		if (evt.type()==IntegratorNonintervalEvent.START){
			openFile();
		}
		else if (evt.type()==IntegratorNonintervalEvent.DONE){
			closeFile();
		}

	}

	private AfterPBC afterPBCinstance;
	private Phase phase;
	//private AtomIteratorLeafAtoms iterator;
	private AtomIteratorMolecule iterator;
	private int writeInterval;
	private int intervalCount;
	private String fileName;
	private FileWriter fileWriter;
	private int counter;
	
	/*
	 * -------------------------SUBCLASS AfterPBC----------------------------------
	 */
	
	private static class AfterPBC implements IntegratorIntervalListener{
		
		public AfterPBC(Space space,Species[] species){
			workVector = space.makeVector();
			iterator = new AtomIteratorMolecule(species);
		}
		
		// Method called in main class (see above)
		public int [][] getAtomPBIarray(){
			return atomPBIarray;
		}
				
		public void setPhase(Phase phase){
			atomOldCoord = new Vector[phase.atomCount()];
			for (int j=0; j < atomOldCoord.length; j++){
				atomOldCoord[j] = phase.space().makeVector();
			}
						
			atomPBIarray = new int[phase.atomCount()][phase.space().D()];
			iterator.setPhase(phase);
			phaseDim = phase.getBoundary().getDimensions();
			updateAtomOldCoord();
		}
		
		// Method called in main and sub class (see directly above and above)
		public void updateAtomOldCoord(){
			iterator.reset();
			int i=0;
			while (iterator.hasNext()){
				atomOldCoord[i].E(((AtomLeaf)iterator.nextAtom()).coord.position());
				i++;
			}
		}
		
		public void intervalAction(IntegratorIntervalEvent evt) {
			iterator.reset();
			int i=0;
			
			// workVector is modified to hold a value of box lengths an atom has traveled
			// atomPBIarray is filled here
			while (iterator.hasNext()){
				workVector.E(atomOldCoord[i]);
				workVector.ME(((AtomLeaf)iterator.nextAtom()).coord.position()); 
				workVector.DE(phaseDim);
				
				for (int j=0;j < phaseDim.D();j++){
					
					// Before Math.round, workVector is -/+ 0.9999,1.000,1.0001,0.000
					// Value will truncate when added to atomPBIarray, we must make workVector a whole number
					atomPBIarray[i][j] += Math.round(workVector.x(j));
					
				}
				i++;
			}
		}
		
		// **
		public int getPriority() {
			return 200;
		}
		
		private Vector phaseDim;
		private int [][] atomPBIarray;
		private Vector workVector;
		private Vector [] atomOldCoord;
		//private AtomIteratorLeafAtoms iterator;
		private AtomIteratorMolecule iterator;
	}
	
}
