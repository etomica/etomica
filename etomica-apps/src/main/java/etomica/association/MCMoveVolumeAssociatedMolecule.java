/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import etomica.action.AtomActionTranslateBy;
import etomica.action.MoleculeChildAtomAction;
import etomica.action.WriteConfiguration;
import etomica.atom.IAtomOriented;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveMolecular;
import etomica.models.OPLS.SpeciesAceticAcid;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.molecule.MoleculeArrayList;
import etomica.molecule.iterator.MoleculeIterator;
import etomica.molecule.iterator.MoleculeIteratorAllMolecules;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Pressure;
import etomica.util.random.IRandom;

import java.io.FileWriter;
import java.io.IOException;

/**
 * Monte Carlo volume-change move for associating fluids simulations in the NPT ensemble.
 *
 * @author Hye Min Kim
 */
public class MCMoveVolumeAssociatedMolecule extends MCMoveBoxStep implements MoleculeAgentSource, MCMoveMolecular {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    private final int D;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected final MoleculeIteratorAllMolecules moleculeIterator;
    protected final Vector r;
    protected final Vector dr, dr2;
    private transient double uOld, hOld, vNew, vScale, hNew;
    private transient double uNew = Double.NaN;
    protected AssociationManagerMolecule associationManager;
    protected int numAssociatedMolecules;
    protected int numMer;
    protected MoleculeAgentManager moleculeAgentManager;
    protected final MoleculeArrayList smerList;
    public static boolean dodebug;
    protected FileWriter fileWriter;
    protected IAssociationHelperMolecule associationHelper;
    protected Vector groupTranslationVector;
    protected MoleculeChildAtomAction moveMoleculeAction;
    protected final Simulation sim;

    public MCMoveVolumeAssociatedMolecule(Simulation sim, PotentialMaster potentialMaster,
                                          Space _space) {
        this(sim, potentialMaster, sim.getRandom(), _space, 1.0);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeAssociatedMolecule(Simulation sim, PotentialMaster potentialMaster, IRandom random,
                                          Space _space, double pressure) {
        super(potentialMaster);
        this.random = random;
        this.sim = sim;
        smerList = new MoleculeArrayList();
        this.D = _space.D();
        r = _space.makeVector();
        this.dr = _space.makeVector();
        this.dr2 = _space.makeVector();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        setStepSizeMax(1.0);
        setStepSizeMin(0.0);
        setStepSize(0.10);
        setPressure(pressure);
        energyMeter.setIncludeLrc(true);
        affectedAtomIterator = new AtomIteratorLeafAtoms();
        AtomActionTranslateBy translator = new AtomActionTranslateBy(_space);
        groupTranslationVector = translator.getTranslationVector();
        moveMoleculeAction = new MoleculeChildAtomAction(translator);
        moleculeIterator = new MoleculeIteratorAllMolecules();
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        affectedAtomIterator.setBox(p);
        moleculeIterator.setBox(p);
        moleculeAgentManager = new MoleculeAgentManager(sim, box, this);
        
    }
    
    public void setAssociationManager(AssociationManagerMolecule associationManager, IAssociationHelperMolecule associationHelper) {
    	this.associationManager = associationManager;
    	this.associationHelper = associationHelper;
    }
    
    public boolean doTrial() {
        double vOld = box.getBoundary().volume();
        uOld = energyMeter.getDataAsScalar();
        if(uOld > 1e8) {
            throw new RuntimeException("atom "+" in box "+box+" has an overlap");
        }
        hOld = uOld + pressure*vOld;
        vScale = (2.*random.nextDouble()-1.)*stepSize;
        vNew = vOld * Math.exp(vScale); //Step in ln(V)
        double rScale = Math.exp(vScale/D);
        numAssociatedMolecules = associationManager.getAssociatedMolecules().getMoleculeCount();
        if (dodebug){
        	WriteConfiguration writeConfiguration = new WriteConfiguration(Space3D.getInstance());
        	writeConfiguration.setBox(box);
        	writeConfiguration.setConfName("old");
        	writeConfiguration.actionPerformed();
        }
        scaleAtoms(rScale);
        if (dodebug){
        	System.out.println("rScale "+rScale);
        	WriteConfiguration writeConfiguration = new WriteConfiguration(Space3D.getInstance());
        	writeConfiguration.setBox(box);
        	writeConfiguration.setConfName("new");
        	writeConfiguration.actionPerformed();
        	throw new RuntimeException();
        }
        uNew = energyMeter.getDataAsScalar();
        hNew = uNew + pressure*vNew;
        return true;
    }//end of doTrial
    
    protected void scaleAtoms(double rScale) {
    	if (dodebug) {
    		try{
    			fileWriter = new FileWriter("smerList");
    		}catch(IOException e){
    			throw new RuntimeException(e);
    		}
    	}
    	numMer = 0;
    	IMoleculeList moleculeList = box.getMoleculeList();// all atoms in the box
    	for (int i=0; i< moleculeList.getMoleculeCount(); i++){
        	IMolecule molecule = moleculeList.getMolecule(i);
        	((Agent)moleculeAgentManager.getAgent(molecule)).nMolecules = 0;
    	}
    	for (int i=0; i< moleculeList.getMoleculeCount(); i++){
        	IMolecule molecule = moleculeList.getMolecule(i);
        	if (((Agent)moleculeAgentManager.getAgent(molecule)).nMolecules!= 0){
        		continue;
        	}
        	numMer++;
        	associationHelper.populateList(smerList,molecule, false);

        	if (dodebug){
        		try{
        			fileWriter.write(smerList.toString()+"\n");
        		}catch(IOException e){
        			throw new RuntimeException(e);
    			}
        	}
        	for (int j=0; j<smerList.getMoleculeCount(); j+=1){
        		Agent jAgent = ((Agent)moleculeAgentManager.getAgent(smerList.getMolecule(j)));
        		jAgent.nMolecules = smerList.getMoleculeCount();
        		jAgent.nextMolecule = j<smerList.getMoleculeCount()-1 ? smerList.getMolecule(j+1) : null;
        	}
        	if (smerList.getMoleculeCount() == 1){
            	groupTranslationVector.Ea1Tv1(rScale-1.0, positionDefinition(molecule));
                moveMoleculeAction.actionPerformed(molecule);
        		continue;
        	}
        	
        	r.E(0.0);
            for (int j = 0; j<smerList.getMoleculeCount(); j+=1){
            	dr.Ev1Mv2(positionDefinition(smerList.getMolecule(j)), r);//dr = distance from the atom to the center of the mass
            	box.getBoundary().nearestImage(dr);
            	groupTranslationVector.Ev1Pv2(r, dr);
            	groupTranslationVector.ME(positionDefinition(smerList.getMolecule(j)));
                moveMoleculeAction.actionPerformed(smerList.getMolecule(j));
            	r.PEa1Tv1(1.0/(j+1), dr);
            }
            r.TE(rScale-1.0);
            groupTranslationVector.E(r);
            for (int j = 0; j<smerList.getMoleculeCount(); j+=1){
                moveMoleculeAction.actionPerformed(smerList.getMolecule(j));
            }
    	}
    	r.Ea1Tv1(rScale, box.getBoundary().getBoxSize());
        box.getBoundary().setBoxSize(r);//scale the boundary
        if (dodebug){
        	try{
        		fileWriter.close();
			}catch(IOException e){
				throw new RuntimeException(e);
			}
        }
    }
    
    protected void unscaleAtoms(double rScale) {
    	IMoleculeList moleculeList = box.getMoleculeList();
        
        for (int i=0; i< moleculeList.getMoleculeCount(); i++){
        	IMolecule molecule = moleculeList.getMolecule(i);
        	Agent iAgent = (Agent)moleculeAgentManager.getAgent(molecule);
        	if (iAgent.nMolecules == 0){
        		continue;
        	}
        	IMoleculeList bondedMolecules = associationManager.getAssociatedMolecules(molecule);
        	if (bondedMolecules.getMoleculeCount() == 0) {
        		groupTranslationVector.Ea1Tv1(rScale-1.0, positionDefinition(molecule));
                moveMoleculeAction.actionPerformed(molecule);//scale the position directly
        		if (dodebug){
        			if (molecule.getIndex() == 19 || molecule.getIndex() == 76||molecule.getIndex() == 127 || molecule.getIndex() == 139||molecule.getIndex() == 196){
        				System.out.println("One of these is monomer:1st "+molecule);
        				IAtomOriented atom19 = (IAtomOriented)box.getLeafList().getAtom(19);
        				IAtomOriented atom139 = (IAtomOriented)box.getLeafList().getAtom(139);
        				System.out.println(" orientation of atom19: "+atom19.getOrientation().getDirection());
    					System.out.println(" orientation of atom139: "+atom139.getOrientation().getDirection());
        			}
        		}
        		continue;
        	}
        	if (iAgent.nMolecules == 1){
        		iAgent.nMolecules = 0;//prevent unscaling again
        		groupTranslationVector.Ea1Tv1(rScale-1.0, positionDefinition(molecule));
                moveMoleculeAction.actionPerformed(molecule);//scale the position directly
        		if (dodebug){
        			if (molecule.getIndex() == 19 || molecule.getIndex() == 76||molecule.getIndex() == 127 || molecule.getIndex() == 139||molecule.getIndex() == 196){
        				System.out.println("One of these is monomer:2nd "+molecule);
        			}
        		}
        		continue;
        	}
        	r.E(0.0);
        	IMolecule jMolecule = molecule;
        	int count = 0;
        	while (jMolecule != null){
            	dr.Ev1Mv2(positionDefinition(jMolecule), r);//dr = distance from the atom to the center of the mass
            	box.getBoundary().nearestImage(dr);
            	groupTranslationVector.Ev1Pv2(r, dr);
            	groupTranslationVector.ME(positionDefinition(jMolecule));
                moveMoleculeAction.actionPerformed(jMolecule);//move molecule2 to the outside of box
            	r.PEa1Tv1(1.0/(count+1), dr);
            	count++;
            	Agent jAgent = (Agent)moleculeAgentManager.getAgent(jMolecule);
            	jMolecule = jAgent.nextMolecule;
            	if (dodebug){
            		if (molecule.getIndex() == 19){
            			System.out.println("jAtom of atom19= "+jMolecule);
            			System.out.println("dr of atom19= "+dr);
            			System.out.println("r of atom19= "+r);
            		}
            	}
            }
            r.TE(rScale-1.0);//cancel the scale
            if (dodebug){
        		if (molecule.getIndex() == 19){
        			System.out.println("r of atom19 after unscaling= "+r);
        		}
        	}
        	jMolecule = molecule;
        	groupTranslationVector.E(r);
            while (jMolecule != null){
                moveMoleculeAction.actionPerformed(jMolecule);
            	Agent jAgent = (Agent)moleculeAgentManager.getAgent(jMolecule);
            	jAgent.nMolecules = 0;//prevent unscaling again
            	jMolecule = jAgent.nextMolecule;
            }
        }
		r.Ea1Tv1(rScale, box.getBoundary().getBoxSize());
        box.getBoundary().setBoxSize(r);//scale the boundary
    }
    

    
    public double getA() {
    	IMoleculeList moleculeList = box.getMoleculeList();
    	for (int i=0; i< moleculeList.getMoleculeCount(); i++){
        	IMolecule molecule = moleculeList.getMolecule(i);
        	if(associationHelper.populateList(smerList, molecule, true)){
        		return 0.0;//reject
        	}
        	
        	for (int j=0; j<smerList.getMoleculeCount(); j+=1){
        		Agent jAgent = ((Agent)moleculeAgentManager.getAgent(molecule));
        		if (jAgent.nMolecules != smerList.getMoleculeCount()){
        			return 0;
        		}
        	}
    	}
        return Math.exp((numMer+1)*vScale);
    }
    
    
    public double getB() {
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */}
    
    public void rejectNotify() {
    	double rScale = Math.exp(vScale/D);
        unscaleAtoms(1/rScale);
        if (dodebug){
        	WriteConfiguration writeConfiguration = new WriteConfiguration(Space3D.getInstance());
        	writeConfiguration.setBox(box);
        	writeConfiguration.setConfName("newnew");
        	writeConfiguration.actionPerformed();
        }
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }
    
    public Vector positionDefinition(IMolecule molecule){
    	return molecule.getChildList().getAtom(SpeciesAceticAcid.indexC).getPosition();
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}


	public Object makeAgent(IMolecule a) {
		return new Agent();
	}

	public void releaseAgent(Object agent, IMolecule molecule) {
		
	}
	public static class Agent {
		public int nMolecules;//number of molecules in the smer
		public IMolecule nextMolecule;
	}

	public Class getMoleculeAgentClass() {
		return Agent.class;
	}

	public MoleculeIterator affectedMolecules(Box box) {
		return moleculeIterator;
	}

}
