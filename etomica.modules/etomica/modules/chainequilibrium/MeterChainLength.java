/*
 * Created on May 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.modules.chainequilibrium;

import java.io.Serializable;

import etomica.atom.Atom;
import etomica.atom.AtomAgentManager;
import etomica.atom.AtomAgentManager.AgentSource;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.meter.Meter;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.phase.Phase;
import etomica.units.Fraction;
import etomica.units.Null;
import etomica.util.NameMaker;

/**
 * @author Matt Moynihan MoleuclarCount returns an array with the number of
 *         atoms In molecules with [1,2,3,4,5,6,7-10,10-13,13-25, <25] atoms
 */
public class MeterChainLength implements Meter, Serializable, AgentSource {

    public MeterChainLength(ReactionEquilibrium sim) {
        setName(NameMaker.makeName(this.getClass()));
        agentSource = sim;
        setupData(40);
    }

    /**
     * Creates the data object (a DataFunction) to be returned by getData().  
     * data wraps the histogram's double[] so copying is not needed.
     */
    protected void setupData(int maxChainLength) {

        DataInfo xDataInfo = new DataInfo("Chain Length Distribution",Fraction.DIMENSION, DataDoubleArray.getFactory(new int[]{maxChainLength}));
        data = new DataFunction(new int[]{maxChainLength});
        dataInfo = new DataInfoFunction("Chain Length Distribution", Null.DIMENSION, new int[]{maxChainLength}, new DataInfo[]{xDataInfo});
        double[] x = data.getXData(0).getData();
        for (int i=0; i<maxChainLength; i++) {
            x[i] = i+1;
        }
    }
    
    public Class getAgentClass() {
        return AtomTag.class;
    }
    
    public Object makeAgent(Atom a) {
        return new AtomTag();
    }
    
    // does nothing
    public void releaseAgent(Object agent, Atom atom) {}

    //returns the number of molecules with [1,2,3,4,5,6,7-10,10-13,13-25, >25]
    // atoms
    public Data getData() {
        agents = agentSource.getAgents(phase);
        atomTags = (AtomTag[])tagManager.getAgents();
        
        double[] histogram = data.getData();
        for (int i=0; i<histogram.length; i++) {
            histogram[i] = 0;
        }

        // untag all the Atoms
        iterator.reset();
        while (iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            atomTags[a.getGlobalIndex()].tagged = false;
        }

        iterator.reset();

        while(iterator.hasNext()) {
            Atom a = iterator.nextAtom();
            // if an Atom is tagged, it was already counted as part of 
            // another chain
            if (atomTags[a.getGlobalIndex()].tagged) continue;

            int chainLength = recursiveTag(a);
            
            if (chainLength-1 < histogram.length) {
                histogram[chainLength-1] += chainLength;
            }
            else {
                histogram[histogram.length-1] += chainLength;
            }
        }

        for (int i=0; i<histogram.length; i++) {
            histogram[i] /= phase.atomCount();
        }
        
        return data;
    }

    protected int recursiveTag(Atom a) {
        atomTags[a.getGlobalIndex()].tagged = true;

        Atom[] nbrs = agents[a.getGlobalIndex()];

        int ctr = 1;
        
        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(atomTags[nbrs[i].getGlobalIndex()].tagged) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            ctr += recursiveTag(nbrs[i]);
        }
        return ctr;
        
    }
    
    public Phase getPhase() {
        return phase;
    }
    
    public void setPhase(Phase phase) {
        this.phase = phase;
        if (tagManager != null) {
            // allow old agentManager to de-register itself as a PhaseListener
            tagManager.setPhase(null);
        }
        tagManager = new AtomAgentManager(this,phase);
        
        iterator.setPhase(phase);
    }

    public DataInfo getDataInfo() {
        return dataInfo;
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    private final AtomIteratorLeafAtoms iterator = new AtomIteratorLeafAtoms();
    private Phase phase;
    private String name;
    private AtomAgentManager tagManager;
    private AtomTag[] atomTags;
    private final ReactionEquilibrium agentSource;
    private Atom[][] agents;
    private DataFunction data;
    private DataInfoFunction dataInfo;
    
    public static class AtomTag {
        public boolean tagged;
    }

}