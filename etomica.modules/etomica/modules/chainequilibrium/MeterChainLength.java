/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.phase.Phase;
import etomica.units.Fraction;
import etomica.util.NameMaker;

/**
 * @author Matt Moynihan MoleuclarCount returns an array with the number of
 *         atoms In molecules with [1,2,3,4,5,6,7-10,10-13,13-25, <25] atoms
 */
public class MeterChainLength implements Meter, Atom.AgentSource {

    private final DataDoubleArray data;
    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private Phase phase;
    private String name;
    private AtomAgentManager tagManager;
    private AtomTag[] atomTags;
    private final ReactionEquilibrium agentSource;
    private Atom[][] agents;

    public MeterChainLength(ReactionEquilibrium sim) {
        data = new DataDoubleArray("Chain Length Distribution",Fraction.DIMENSION,10);
        setName(NameMaker.makeName(this.getClass()));
        agentSource = sim;
    }

    public Object makeAgent(Atom a) {
        //System.out.println("added atom Tag");
        return new AtomTag();
    }
    
    public void releaseAgent(Object agent) {}

    //returns the number of molecules with [1,2,3,4,5,6,7-10,10-13,13-25, >25]
    // atoms
    public Data getData() {
        agents = agentSource.getAgents(phase);
        atomTags = (AtomTag[])tagManager.getAgents();
        
        //System.out.println("Meter Dimer Frac: ran MolculeFraction Method");
        //double[] array = new double[1];
        //array[0] = 2.1;
        data.E(0);
        
        iterator.reset();
        while (iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            atomTags[a.getGlobalIndex()].reset();
        }

        iterator.reset();
        double[] x = data.getData();

        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            if (!atomTags[a.getGlobalIndex()].white()) continue;

            int j = recursive(a);
            if ((j < 7) || (j == 7)) {
                x[(j - 1)] += j;
            } else {
                if ((j == 8) || (j == 9) || (j == 10)) {
                    x[6] += j;
                }
                if ((j == 11) || (j == 12) || (j == 13)) {
                    x[7] += j;
                } else {
                    if (j < 25) {
                        x[8] += j;
                    } else {
                        x[9] += j;
                    }
                }
            }
        }
        data.TE(1.0/phase.atomCount());
        return data;
    }

    //this function tells you if next is terminal (I.e it is bonded once)
    public boolean terminal(Atom current, Atom next) {
        int j = agents[next.getGlobalIndex()].length; //check INDEXING
        for (int i = 0; i != j; ++i) {
            if (agents[next.getGlobalIndex()][i] == current) {

            } else {
                if (agents[next.getGlobalIndex()][i] != null) {
                    return false;
                }
            }
        }
        return true;
    }

    //returns true only if there is NOTHING bonded to the atom
    public boolean empty(Atom a) {
        int j = agents[a.getGlobalIndex()].length; //check INDEXING
        for (int i = 0; i != j; ++i) {
            if (agents[a.getGlobalIndex()][i] != null) {
                return false;
            }
        }
        return true;
    }

    //returns the number of unmarked atoms it is bonded too (not including
    // itself)
    public int NumberUmarkedBonded(Atom a) {
        int j = agents[a.getGlobalIndex()].length; //check INDEXING
        int ctr = 0;

        // Loops through all the spots atom a can use to bond
        for (int i = 0; i != j; ++i) {
            Atom current = agents[a.getGlobalIndex()][i];
            if (current != null) { // There is an atom
                if (atomTags[current.getGlobalIndex()].white()) {
                    ctr++;
                }
            }
        }
        return ctr;
    }

    public int NumberBonded(Atom a) {
        //System.out.println("Molecular Count: ran Number Bonded");
        int j = agents[a.getGlobalIndex()].length;
        int ctr = 0;
        for (int i = 0; i != j; ++i) {
            Atom current = agents[a.getGlobalIndex()][i];
            if (current != null) { // There is an atom
                if (current != a) {
                    ctr++;
                }
            }
        }
        return ctr;
    }

    //This returns an array of all the white atoms attached to atom a
    public Atom[] findAllWhite(Atom a) {
        //System.out.println("Molecular Count: find all white method");
        int j = agents[a.getGlobalIndex()].length; //check INDEXING
        Atom[] array = new Atom[j];
        int ctr = 0;

        // Sets all spots in the new array equal to NULL
        for (int i = 0; i != j; ++i) {
            array[i] = null;
        }

        // Loops through all the spots atom a can use to bond
        for (int i = 0; i != j; ++i) {
            Atom current = agents[a.getGlobalIndex()][i];
            if (current != null) { // There is an atom
                if (atomTags[current.getGlobalIndex()].white()) {
                    array[ctr] = current;
                    ++ctr;
                }
            }
        }
        return array;
    }

    public int recursive(Atom a) {
        if (atomTags[a.getGlobalIndex()].white()) {
            //System.out.println("Ran recursive");
            atomTags[a.getGlobalIndex()].process();
        }
        Atom[] nbrs = agents[a.getGlobalIndex()];

        int ctr = 1;
        
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(!atomTags[nbrs[i].getGlobalIndex()].white()) continue; 
            ctr += recursive(nbrs[i]);
        }
        return ctr;
        //int j = ((Atom[]) a.allatomAgents[nbrIndex]).length;
        
    }
    
    
    /* (non-Javadoc)
     * @see etomica.data.meter.Meter#getPhase()
     */
    public Phase getPhase() {
        return phase;
    }
    
    /* (non-Javadoc)
     * @see etomica.data.meter.Meter#setPhase(etomica.phase.Phase)
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
        if (tagManager != null) {
            // allow old agentManager to de-register itself as a PhaseListener
            tagManager.setPhase(null);
        }
        tagManager = new AtomAgentManager(this,phase);
        
        iterator.setPhase(phase);
    }

    /* (non-Javadoc)
     * @see etomica.data.DataSource#getDataInfo()
     */
    public DataInfo getDataInfo() {
        return data.getDataInfo();
    }
    
    /* (non-Javadoc)
     * @see etomica.EtomicaElement#getName()
     */
    public String getName() {
        return name;
    }
    /* (non-Javadoc)
     * @see etomica.EtomicaElement#setName(java.lang.String)
     */
    public void setName(String name) {
        this.name = name;
    }
}