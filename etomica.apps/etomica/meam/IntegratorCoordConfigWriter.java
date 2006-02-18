package etomica.meam;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.integrator.Integrator;
import etomica.integrator.IntegratorIntervalEvent;
import etomica.integrator.IntegratorIntervalListener;
import etomica.integrator.IntegratorNonintervalEvent;
import etomica.integrator.IntegratorNonintervalListener;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.space.Vector;

public class IntegratorCoordConfigWriter implements IntegratorIntervalListener,
		IntegratorNonintervalListener {
	
	
	public IntegratorCoordConfigWriter(Space space, String fileName){
		
		afterPBCinstance = new AfterPBC(space);
		
		this.fileName = fileName;
		iterator = new AtomIteratorLeafAtoms();
		setWriteInterval(1);
	}
	
	public void setPhase(Phase phase){
		
		this.phase = phase;
		iterator.setPhase(phase);
		afterPBCinstance.setPhase(phase);
	}
	
	public void setIntegrator(Integrator integrator){
		
		integrator.addListener(this);
		integrator.addListener(afterPBCinstance);
	}
	
	public void openFile(){
		try { 
			fileWriter = new FileWriter(fileName, false);
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
	
	public void setWriteInterval(int writeInterval){
		this.writeInterval = writeInterval;
		intervalCount = writeInterval;
	}
	
	public void intervalAction(IntegratorIntervalEvent evt) {
		afterPBCinstance.updateAtomOldCoord();
		if (--intervalCount == 0){
			Vector phasedim = phase.getBoundary().getDimensions();
			int [][] atomPBIarray = afterPBCinstance.getAtomPBIarray();
					
			try {
				iterator.reset();
				int i=0;
				
				while (iterator.hasNext()){
					Vector atomPosition = ((AtomLeaf)iterator.nextAtom()).coord.position();
					for (int j=0;j < phasedim.D();j++){
						double actualDistance;
						actualDistance = atomPBIarray[i][j] * phasedim.x(j) + atomPosition.x(j);
						fileWriter.write("     "+actualDistance);
					}
				fileWriter.write("\n");
				i++;
				}
			}
			catch (IOException e) {
	            throw new RuntimeException(e);
	        }
			intervalCount = writeInterval;
		}
	}
		

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
	private AtomIteratorLeafAtoms iterator;
	private int writeInterval;
	private int intervalCount;
	private String fileName;
	private FileWriter fileWriter;
	
	private static class AfterPBC implements IntegratorIntervalListener{
		
		public AfterPBC(Space space){
			workVector = space.makeVector();
			iterator = new AtomIteratorLeafAtoms();
		}
		
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
			while (iterator.hasNext()){
				workVector.E(((AtomLeaf)iterator.nextAtom()).coord.position());
				workVector.ME(atomOldCoord[i]);
				workVector.DE(phaseDim);
				for (int j=0;j < phaseDim.D();j++){
					atomPBIarray[i][j] += workVector.x(j);
				}
				i++;
			}
		}

		public int getPriority() {
			return 200;
		}
		
		private Vector phaseDim;
		private int [][] atomPBIarray;
		private Vector workVector;
		private Vector [] atomOldCoord;
		private AtomIteratorLeafAtoms iterator;
		
	}
	
}
