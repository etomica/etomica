package etomica.action;

import java.io.FileWriter;
import java.io.IOException;

import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;
import etomica.space.IVector;

/**
 * Dumps a phase's configuration to a file.  The coordinates are written in a 
 * format that can be read in by ConfigurationFile.  The output file has a 
 * "pos_new" extension, which should be renamed to "pos" for use with
 * ConfigurationFile.
 */
public class WriteConfiguration implements Action {

    /**
     * Returns a default instance.  The configuration name and phase should be 
     * set via setConfName and setPhase before use.
     */
    public WriteConfiguration() {
        iterator = new AtomIteratorLeafAtoms(null);
    }

    /**
     * Sets the configuration name.  The file written to is newConfName.pos_new
     */
    public void setConfName(String newConfName) {
        confName = newConfName;
    }
    
    /**
     * Returns the configuration name.  The file written to is confName.pos_new
     */
    public String getConfName() {
        return confName;
    }
    
    /**
     * Sets the phase whose atom coordinates get written to the file.
     */
    public void setPhase(Phase newPhase) {
        phase = newPhase;
        iterator.setPhase(newPhase);
        setDoApplyPBC(true);
    }
    
    /**
     * Returns the phase whose atom coordinates get written to the file.
     */
    public Phase getPhase() {
        return phase;
    }
    
    /**
     * Directs the writer to apply periodic boundary conditions or not (true 
     * by default).
     */
    public void setDoApplyPBC(boolean newDoApplyPBC) {
        doApplyPBC = newDoApplyPBC;
    }
    
    /**
     * Returns true if PBC are applied to coordinates written to the file.
     */
    public boolean getDoApplyPBC() {
        return doApplyPBC;
    }
    
    /**
     * Writes the leaf Atom coordinates to the file confName.pos_new.  If the
     * file exists, it is overwritten.
     */
    public void actionPerformed() {
        FileWriter fileWriter;
        String fileName = confName + ".pos_new";
        try { 
            fileWriter = new FileWriter(fileName);
        }catch(IOException e) {
            System.err.println("Cannot open "+fileName+", caught IOException: " + e.getMessage());
            return;
        }
        try {
            iterator.reset();
            IVector writePosition = phase.getSpace().makeVector();
            while (iterator.hasNext()) {
                AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
                writePosition.E(atom.getPosition());
                if (doApplyPBC) {
                    IVector shift = phase.getBoundary().centralImage(writePosition);
                    if (!shift.isZero()) {
                        writePosition.PE(shift);
                    }
                }
                
                fileWriter.write(writePosition.x(0)+"");
                for (int i=1; i<writePosition.getD(); i++) {
                	fileWriter.write(" "+writePosition.x(i));
                }
                fileWriter.write("\n");
            }
            fileWriter.close();
        } catch(IOException e) {
            System.err.println("Problem writing to "+fileName+", caught IOException: " + e.getMessage());
        }
    }

    private String confName;
    private Phase phase;
    private boolean doApplyPBC;
    private final AtomIteratorLeafAtoms iterator;
}
