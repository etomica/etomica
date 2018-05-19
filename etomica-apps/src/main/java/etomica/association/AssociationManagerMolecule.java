/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.atom.AtomArrayList;
import etomica.box.Box;
import etomica.box.BoxAgentManager;
import etomica.integrator.mcmove.MCMoveEvent;
import etomica.integrator.mcmove.MCMoveMolecular;
import etomica.integrator.mcmove.MCMoveTrialCompletedEvent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.iterator.MoleculeIterator;
import etomica.nbr.cell.molecule.Mpi1ACell;
import etomica.nbr.cell.molecule.NeighborCellManagerMolecular;
import etomica.simulation.Simulation;
import etomica.util.IEvent;
import etomica.util.IListener;

/**
 * Class to define and track molecule associations.  Constructed given an iterator
 * that defines the set of molecules that is managed, and an association definition
 * that is used to determine if molecules are associated.  
 * To incorporate an instance of this class in a simulation, add an instance of 
 * this EnergySum inner class to the simulation:<br>
 *       sim.setEnergySum(associationManager.new EnergySum());
 *
 * If a MC simulation, this should be done before any MCMove classes are instantiated.
 * Also, register this as a listener to IntegratorMC:<br>
 *       integrator.addMCMoveListener(associationManager);
 *       
 *@author Hye Min Kim
 */
public class AssociationManagerMolecule implements MoleculeAgentSource,IListener {
    
    private AssociationDefinitionMolecule associationDefinition;
    private final Box box;
    private final MoleculeAgentManager agentManager;
    private final Mpi1ACell neighborIterator;
    private final MoleculeArrayList associatedMolecules = new MoleculeArrayList();
    private final IListener mcMoveListener;

    public AssociationManagerMolecule(Simulation sim, Box box, BoxAgentManager<NeighborCellManagerMolecular> cellAgentManager, AssociationDefinitionMolecule definition, double range) {
    	this.box = box;
    	agentManager = new MoleculeAgentManager(sim,box,this);
        associationDefinition = definition;
        this.neighborIterator = new Mpi1ACell(3,range,cellAgentManager);
        mcMoveListener = cellAgentManager.getAgent(box).makeMCMoveListener(); 
        cellAgentManager.getAgent(box).setDoApplyPBC(true);
    }
    
    public AssociationDefinitionMolecule getAssociationDefinition() {
    	return associationDefinition;
    }
    public void initialize() {
        IMoleculeList moleculeList = box.getMoleculeList();//list of all atoms in this box
        for (int i = 0; i<moleculeList.size(); i+=1) {
        	IMolecule moleculei = moleculeList.get(i);
        	((MoleculeArrayList)agentManager.getAgent(moleculei)).clear();
        }
        for (int i = 0; i<moleculeList.size()-1; i+=1) {
        	IMolecule moleculei = moleculeList.get(i);//definition of atom i
        	for (int j = i+1; j<moleculeList.size(); j+=1) {
            	IMolecule moleculej = moleculeList.get(j);
            	if(associationDefinition.isAssociated(moleculei,moleculej)) {
                    ((MoleculeArrayList)agentManager.getAgent(moleculei)).add(moleculej);//i and j are associated
                    ((MoleculeArrayList)agentManager.getAgent(moleculej)).add(moleculei);
                }
        	}
         }
        associatedMolecules.clear();
        for (int i = 0; i<moleculeList.size(); i+=1) {
        	IMolecule moleculei = moleculeList.get(i);
        	if(((MoleculeArrayList)agentManager.getAgent(moleculei)).size() > 0){
        		associatedMolecules.add(moleculei);
        	}
        }
                
    }
    
    public IMoleculeList getAssociatedMolecules() {return associatedMolecules;}
    
    //need also to handle associatedAtoms list
    public void actionPerformed(IEvent evt) {
    	mcMoveListener.actionPerformed(evt);
    	MCMoveEvent mcEvent = (MCMoveEvent)evt;
        if(mcEvent instanceof MCMoveTrialCompletedEvent && ((MCMoveTrialCompletedEvent) mcEvent).isAccepted()) {
        	return;
            }
        
        MoleculeIterator iterator = ((MCMoveMolecular)mcEvent.getMCMove()).affectedMolecules(box);
        iterator.reset();
        neighborIterator.setBox(box);
        for (IMolecule moleculei = iterator.nextMolecule();moleculei != null; moleculei =iterator.nextMolecule()){
        	MoleculeArrayList listi = (MoleculeArrayList)agentManager.getAgent(moleculei);//list of the old bonds
        	if (listi.size() > 0){
        		for (int i = 0; i<listi.size(); i++){
            		IMolecule moleculej = listi.get(i);
            		MoleculeArrayList listj = (MoleculeArrayList)agentManager.getAgent(moleculej);
        			listj.remove(listj.indexOf(moleculei));//remove atom i from the listj
        			if ( listj.size() == 0) {
        				associatedMolecules.remove(associatedMolecules.indexOf(moleculej));
        			}
            	}
            	listi.clear();//remove all the elements from listi
            	associatedMolecules.remove(associatedMolecules.indexOf(moleculei));
        	}
        	
        	neighborIterator.setTarget(moleculei);
        	neighborIterator.reset();
        	for (IMoleculeList moleculeij = neighborIterator.next();moleculeij != null; moleculeij =neighborIterator.next()){
            	IMolecule moleculej = moleculeij.get(0);
            	if (moleculej == moleculei){
            		moleculej = moleculeij.get(1);
            	}
            	if (moleculeij == null) throw new RuntimeException();
            	if (associationDefinition.isAssociated(moleculei, moleculej)){ //they are associated
        			if (listi.size() == 0) {
        				associatedMolecules.add(moleculei);
        			}
        			listi.add(moleculej); //make atom i and atom j to be associated
        			MoleculeArrayList listj = (MoleculeArrayList)agentManager.getAgent(moleculej);
        			if (listj.size() == 0) {
        				associatedMolecules.add(moleculej);
        			}
        			listj.add(moleculei);//make atom i and atom j to be associated

            	}
        	}
        }
        for (IMolecule moleculei = iterator.nextMolecule();moleculei != null; moleculei =iterator.nextMolecule()){
        if (getAssociatedMolecules(moleculei).size() > 1) {
        	//System.out.println("moleculei" +moleculei);
        } else if (getAssociatedMolecules(moleculei).size() == 1){
        	IMolecule moleculej = getAssociatedMolecules(moleculei).get(0);
        	if (getAssociatedMolecules(moleculej).size() == 0){
        		System.out.println("Wrong");
        	} else if(getAssociatedMolecules(moleculej).size() > 1){
        		AtomArrayList listj = (AtomArrayList)agentManager.getAgent(moleculej);
        		System.out.println("Wrong:smer");
        		System.out.println("listj = "+listj+" moleculei= " +moleculei+" moleculej= "+moleculej);
        	} 
        }
        }
    }
    
        
    /**
     * Returns the number of atoms on the list of associations of the given atom.
     */
    public IMoleculeList getAssociatedMolecules(IMolecule molecule) {
       return (MoleculeArrayList)agentManager.getAgent(molecule);
    }

	public Object makeAgent(IMolecule a) {
		return new MoleculeArrayList();
	}

	public void releaseAgent(Object agent, IMolecule atom) {
		
	}
}
