/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.association;

import java.io.FileWriter;
import java.io.IOException;

import etomica.action.WriteConfiguration;
import etomica.api.*;
import etomica.atom.*;
import etomica.box.Box;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.Dimension;
import etomica.units.Pressure;

/**
 * Standard Monte Carlo volume-change move for simulations in the NPT ensemble.
 *
 * @author David Kofke
 */
public class MCMoveVolumeAssociated extends MCMoveBoxStep implements AtomLeafAgentManager.AgentSource<MCMoveVolumeAssociated.Agent> {
    
    private static final long serialVersionUID = 2L;
    protected double pressure;
    private MeterPotentialEnergy energyMeter;
    private final int D;
    private IRandom random;
    protected final AtomIteratorLeafAtoms affectedAtomIterator;
    protected final Vector r;
    protected final Vector dr, dr2;
    private transient double uOld, hOld, vNew, vScale, hNew;
    private transient double uNew = Double.NaN;
    protected AssociationManager associationManager;
    protected int numAssociatedAtoms;
    protected int numMer;//number of molecules
    protected AtomLeafAgentManager<Agent> atomLeafAgentManager;
    protected final AtomArrayList smerList;
    public static boolean dodebug;
    protected FileWriter fileWriter;

    public MCMoveVolumeAssociated(Simulation sim, PotentialMaster potentialMaster,
                                  Space _space) {
        this(potentialMaster, sim.getRandom(), _space, 1.0);
    }
    
    /**
     * @param potentialMaster an appropriate PotentialMaster instance for calculating energies
     * @param space the governing space for the simulation
     */
    public MCMoveVolumeAssociated(PotentialMaster potentialMaster, IRandom random,
                                  Space _space, double pressure) {
        super(potentialMaster);
        this.random = random;
        smerList = new AtomArrayList();
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
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
        affectedAtomIterator.setBox(p);
        atomLeafAgentManager = new AtomLeafAgentManager<Agent>(this,box,Agent.class);
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
        numAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
//        System.out.println("doTrial, associated atoms = " +associationManager.getAssociatedAtoms());//list of the associated atoms
//        associationManager.initialize();
//        if (numAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){//numAssociatedAtoms after move
//        	System.out.println("doTrial, associated atoms .....= " +associationManager.getAssociatedAtoms());//list of the associated atoms
//    		throw new RuntimeException("Oops");
//    	}
//        if (numAssociatedAtoms %2 == 1){
//        	IAtomList list = associationManager.getAssociatedAtoms();
//        	for ( int i = 0; i< numAssociatedAtoms;i++){
//        		System.out.println("list "+list.getAtom(i)+":"+ associationManager.getAssociatedAtoms(list.getAtom(i)));
//        	}
//        	list = box.getLeafList();
//        	IAtomOriented a = (IAtomOriented)list.getAtom(354);
//        	System.out.println("atom 354 " +a.getPosition()+ " " +a.getOrientation().getDirection());
//        	a = (IAtomOriented)list.getAtom(68);
//        	System.out.println("atom 68 " +a.getPosition()+ " " +a.getOrientation().getDirection());
//        	a = (IAtomOriented)list.getAtom(409);
//        	System.out.println("atom 409 " +a.getPosition()+ " " +a.getOrientation().getDirection());
//        	System.out.println("bonded" +associationManager.getAssociationDefinition().isAssociated(list.getAtom(354), list.getAtom(68)));
//        	System.out.println("bonded" +associationManager.getAssociationDefinition().isAssociated(list.getAtom(354), list.getAtom(409)));
//        	System.out.println("bonded" +associationManager.getAssociationDefinition().isAssociated(list.getAtom(68), list.getAtom(409)));
//        	P2HardAssociationCone p = new P2HardAssociationCone(Space3D.getInstance(), 1.0, 1.0, 6.0, 16.0);
//        	p.setBox(box);
//        	System.out.println("energy 354-68 = "+p.energy(new AtomPair(list.getAtom(354), list.getAtom(68))));
//        	System.out.println("energy 354-409 = "+p.energy(new AtomPair(list.getAtom(354), list.getAtom(409))));
//        	System.out.println("energy 68-409 = "+p.energy(new AtomPair(list.getAtom(68), list.getAtom(409))));
//        	throw new RuntimeException("***");
//        }
        //System.out.println("rScale = "+rScale);
        if (dodebug){
        	WriteConfiguration writeConfiguration = new WriteConfiguration(Space3D.getInstance());
        	writeConfiguration.setBox(box);
        	writeConfiguration.setConfName("old");
        	writeConfiguration.actionPerformed();
        }
        scaleAtoms(rScale);//call the method
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
    	IAtomList atomList = box.getLeafList();// all atoms in the box
    	for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	atomLeafAgentManager.getAgent(atom).nAtoms = 0;
    	}
    	for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	if (atomLeafAgentManager.getAgent(atom).nAtoms!= 0){
        		continue;
        	}
        	numMer++;
        	if (populateList(smerList,atom) == 0){
        		System.out.println("smerList: "+smerList);
        		throw new RuntimeException();
        	}
        	if (dodebug){
        		try{
        			fileWriter.write(smerList.toString()+"\n");
        		}catch(IOException e){
        			throw new RuntimeException(e);
    			}
        	}
        	for (int j=0; j<smerList.getAtomCount(); j+=1){
        		Agent jAgent = atomLeafAgentManager.getAgent(smerList.getAtom(j));
        		jAgent.nAtoms = smerList.getAtomCount();
        		jAgent.nextAtom = j<smerList.getAtomCount()-1 ? smerList.getAtom(j+1) : null;//if condition is true, 1st value, if condition is false, 2nd value
        	}
        	//System.out.println("atom= "+atom+" bonded atoms= "+bondedAtoms);
        	if (smerList.getAtomCount() == 1){
//        		if (atom.getLeafIndex() == 226 || atom.getLeafIndex() == 66) {
//        			System.out.println("monomer atom1 before = "+ atom+ "its position = "+atom.getPosition());
//        		}
        		atom.getPosition().TE(rScale);
//        		if (atom.getLeafIndex() == 226 || atom.getLeafIndex() == 66) {
//        			System.out.println("monomer atom1 after = "+ atom+ "its position = "+atom.getPosition());
//        		}
        		continue;
        	}
        	
        	r.E(0.0);
            for (int j = 0; j<smerList.getAtomCount(); j+=1){
//            	if (smerList.getAtom(j).getLeafIndex() == 226 || smerList.getAtom(j).getLeafIndex() == 66) {
//        			System.out.println("atom1 = "+ smerList.getAtom(j)+ "its position = "+smerList.getAtom(j).getPosition());
//        		}
            	dr.Ev1Mv2((smerList.getAtom(j)).getPosition(), r);//dr = distance from the atom to the center of the mass
            	box.getBoundary().nearestImage(dr);
            	smerList.getAtom(j).getPosition().Ev1Pv2(r, dr);//move atom2 to the outside of box
            	r.PEa1Tv1(1.0/(j+1), dr);
            }
            r.TE(rScale-1.0);
            for (int j = 0; j<smerList.getAtomCount(); j+=1){
            	smerList.getAtom(j).getPosition().PE(r);
//            	if (smerList.getAtom(j).getLeafIndex() == 226 || smerList.getAtom(j).getLeafIndex() == 66) {
//            		System.out.println("rScale= "+rScale+"translation scaling "+ r);
//            	}
            }
    	}
    	r.Ea1Tv1(rScale, box.getBoundary().getBoxSize());
    	//System.out.println("box size = " +box.getBoundary().getDimensions());
    	//System.out.println("r = " +r);  	
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
    	IAtomList atomList = box.getLeafList();
        
        for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	Agent iAgent = atomLeafAgentManager.getAgent(atom);
        	if (iAgent.nAtoms == 0){
        		continue;
        	}
        	IAtomList bondedAtoms = associationManager.getAssociatedAtoms(atom);
        	if (bondedAtoms.getAtomCount() == 0) {
//        		if (atom.getLeafIndex() == 226 || atom.getLeafIndex() == 66) {
//        			System.out.println("unscale atom  = "+ atom +" " +atom.getPosition());
//        		}
        		atom.getPosition().TE(rScale);
        		if (dodebug){
        			if (atom.getLeafIndex() == 19 || atom.getLeafIndex() == 76||atom.getLeafIndex() == 127 || atom.getLeafIndex() == 139||atom.getLeafIndex() == 196){
        				System.out.println("One of these is monomer:1st "+atom);
        				IAtomOriented atom19 = (IAtomOriented)box.getLeafList().getAtom(19);
        				IAtomOriented atom139 = (IAtomOriented)box.getLeafList().getAtom(139);
        				System.out.println(" orientation of atom19: "+atom19.getOrientation().getDirection());
    					System.out.println(" orientation of atom139: "+atom139.getOrientation().getDirection());
        			}
        		}
//        		if (atom.getLeafIndex() == 226 || atom.getLeafIndex() == 66) {
//        			System.out.println("unscale atom ... = "+ atom +" " +atom.getPosition());
//        		}
        		continue;//go to another atom
        	}
        	if (iAgent.nAtoms == 1){
        		iAgent.nAtoms = 0;//prevent unscaling again
        		atom.getPosition().TE(rScale);//scale the position directly
        		if (dodebug){
        			if (atom.getLeafIndex() == 19 || atom.getLeafIndex() == 76||atom.getLeafIndex() == 127 || atom.getLeafIndex() == 139||atom.getLeafIndex() == 196){
        				System.out.println("One of these is monomer:2nd "+atom);
        			}
        		}
        		continue;
        	}
        	r.E(0.0);
        	IAtom jAtom = atom;
        	int count = 0;
        	while (jAtom != null){
            	dr.Ev1Mv2(jAtom.getPosition(), r);//dr = distance from the atom to the center of the mass
            	box.getBoundary().nearestImage(dr);
            	jAtom.getPosition().Ev1Pv2(r, dr);//move atom2 to the outside of box
            	r.PEa1Tv1(1.0/(count+1), dr);
            	count++;
            	Agent jAgent = atomLeafAgentManager.getAgent(jAtom);
            	jAtom = jAgent.nextAtom;
            	if (dodebug){
            		if (atom.getLeafIndex() == 19){
            			System.out.println("jAtom of atom19= "+jAtom);
            			System.out.println("dr of atom19= "+dr);
            			System.out.println("r of atom19= "+r);
            		}
            	}
            }
            r.TE(rScale-1.0);//cancel the scale
            if (dodebug){
        		if (atom.getLeafIndex() == 19){
        			System.out.println("r of atom19 after unscaling= "+r);
        		}
        	}
        	jAtom = atom;
            while (jAtom != null){
            	jAtom.getPosition().PE(r);
//            	if (jAtom.getLeafIndex()==226 || jAtom.getLeafIndex()==66){
//            		System.out.println("jAtom= "+jAtom+ "rScale= "+rScale+"translation unscaling "+ r);
//            	}
            	Agent jAgent = atomLeafAgentManager.getAgent(jAtom);
            	jAgent.nAtoms = 0;//prevent unscaling again
            	jAtom = jAgent.nextAtom;
            }
        }
		r.Ea1Tv1(rScale, box.getBoundary().getBoxSize());
        box.getBoundary().setBoxSize(r);//scale the boundary
    }
    protected int populateList(AtomArrayList mySmerList, IAtom atom){
    	mySmerList.clear();
    	mySmerList.add(atom);
    	IAtomList bondList = associationManager.getAssociatedAtoms(atom);
    	double innerRadius = 0.8;
    	double minDistance = 2*(innerRadius*innerRadius)*(1+Math.cos(etomica.units.Degree.UNIT.toSim(27.0)));
    	if (bondList.getAtomCount() > 2){
    		//System.out.println("atom "+ atom+" bondList: "+ bondList);
    		return 0;
    	}
    	if (bondList.getAtomCount() == 2){
    		IAtom atom0 = bondList.getAtom(0);
    		IAtom atom1 = bondList.getAtom(1);
    		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr2);
        	if (dr2.squared() < minDistance){
        		return 0;
        	}
    	}
    	if (bondList.getAtomCount() == 0){
    		return 1;
    	}
    	IAtom thisAtom = bondList.getAtom(0);
    	mySmerList.add(thisAtom);
    	IAtomList bondList1 = associationManager.getAssociatedAtoms(thisAtom);
    	if (bondList1.getAtomCount() > 2){
    		//System.out.println("this atom "+ thisAtom+" bondList1: "+ bondList1);
    		return 0;
    	}
    	if (bondList1.getAtomCount() == 2){
    		IAtom atom0 = bondList1.getAtom(0);
    		IAtom atom1 = bondList1.getAtom(1);
    		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
        	box.getBoundary().nearestImage(dr2);
        	if (dr2.squared() < minDistance){
        		return 0;
        	}
    	}
    	IAtom previousAtom = atom;
    	while (bondList1.getAtomCount() > 1){
    		IAtom nextAtom = bondList1.getAtom(0);
    		if (nextAtom == previousAtom){
    			nextAtom = bondList1.getAtom(1);
    		} 
    		if (nextAtom == atom){
    			return 1;
    		}
    		mySmerList.add(nextAtom);
    		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
    		if (bondList1.getAtomCount() > 2){
    			//System.out.println("next atom "+ nextAtom+" bondList1: "+ bondList1);
        		return 0;
        	}
    		if (bondList1.getAtomCount() == 2){
        		IAtom atom0 = bondList1.getAtom(0);
        		IAtom atom1 = bondList1.getAtom(1);
        		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr2);
            	if (dr2.squared() < minDistance){
            		return 0;
            	}
        	}
    		previousAtom = thisAtom;
    		thisAtom = nextAtom;
    	}
    	if (bondList.getAtomCount()>1){
    		thisAtom = bondList.getAtom(1);
        	mySmerList.add(thisAtom);
        	bondList1 = associationManager.getAssociatedAtoms(thisAtom);
        	if (bondList1.getAtomCount() > 2){
        		//System.out.println("this atom "+ thisAtom+" bondList1: "+ bondList1);
        		return 0;
        	}
        	if (bondList1.getAtomCount() == 2){
        		IAtom atom0 = bondList1.getAtom(0);
        		IAtom atom1 = bondList1.getAtom(1);
        		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
            	box.getBoundary().nearestImage(dr2);
            	
            	if (dr2.squared() < minDistance){
            		return 0;
            	}
        	}
        	previousAtom = atom;
        	while (bondList1.getAtomCount() > 1){
        		IAtom nextAtom = bondList1.getAtom(0);
        		if (nextAtom == previousAtom){
        			nextAtom = bondList1.getAtom(1);
        		} 
        		mySmerList.add(nextAtom);
        		bondList1 = associationManager.getAssociatedAtoms(nextAtom);
        		if (bondList1.getAtomCount() > 2){
        			//System.out.println("next atom "+ nextAtom+" bondList1: "+ bondList1);
            		return 0;
            	}
        		if (bondList1.getAtomCount() == 2){
            		IAtom atom0 = bondList1.getAtom(0);
            		IAtom atom1 = bondList1.getAtom(1);
            		dr2.Ev1Mv2((atom0).getPosition(), (atom1).getPosition());//dr = distance from the atom0 to atom1
                	box.getBoundary().nearestImage(dr2);
                	if (dr2.squared() < minDistance){
                		return 0;
                	}
            	}
        		previousAtom = thisAtom;
        		thisAtom = nextAtom;
        	}
    	}
    	return 1;
    }
    public void setAssociationManager(AssociationManager associationManager) {
    	this.associationManager = associationManager;
    }
    
    public double getA() {
//    	int newnumAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
//    	associationManager.initialize();
//    	if (newnumAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){
//    		throw new RuntimeException("wrong");
//    	}
//    	System.out.println("getA, number of associated atoms = " +associationManager.getAssociatedAtoms().getAtomCount());
//    	System.out.println("number of num.associated atoms = " +numAssociatedAtoms);
//    	if (numAssociatedAtoms > associationManager.getAssociatedAtoms().getAtomCount()){
//    		System.out.println("getA, associated atoms = " +associationManager.getAssociatedAtoms());//list of the associated atoms
//    		throw new RuntimeException ("!!!");
//    	}
    	IAtomList atomList = box.getLeafList();
    	for (int i=0; i< atomList.getAtomCount(); i++){
        	IAtom atom = atomList.getAtom(i);
        	populateList(smerList, atom);
        	for (int j=0; j<smerList.getAtomCount(); j+=1){
        		Agent jAgent = atomLeafAgentManager.getAgent(atom);
        		if (jAgent.nAtoms != smerList.getAtomCount()){
        			return 0;
        		}
        	}
    	}
        return Math.exp((numMer+1)*vScale);
    }
    
    
    public double getB() {
        return -(hNew - hOld);
    }
    
    public void acceptNotify() {  /* do nothing */
    	//System.out.println("accepted volume moving ");
    	//System.out.println("position 388 "+((IAtomPositioned)box.getLeafList().getAtom(388)).getPosition()+"position 115 "+((IAtomPositioned)box.getLeafList().getAtom(115)).getPosition());
    }
    
    public void rejectNotify() {
    	double rScale = Math.exp(vScale/D);
        unscaleAtoms(1/rScale);
        if (dodebug){
        	WriteConfiguration writeConfiguration = new WriteConfiguration(Space3D.getInstance());
        	writeConfiguration.setBox(box);
        	writeConfiguration.setConfName("newnew");
        	writeConfiguration.actionPerformed();
        }
    	//System.out.println("rejected volume moving ");
    	//System.out.println("position 388 "+((IAtomPositioned)box.getLeafList().getAtom(388)).getPosition()+"position 115 "+((IAtomPositioned)box.getLeafList().getAtom(115)).getPosition());
    
//        int newnumAssociatedAtoms = associationManager.getAssociatedAtoms().getAtomCount();
//    	associationManager.initialize();
//    	if (newnumAssociatedAtoms != associationManager.getAssociatedAtoms().getAtomCount()){
//    		System.out.println("newnumAssociatedAtoms= "+newnumAssociatedAtoms);
//    		System.out.println("associationManager.getAssociatedAtoms().getAtomCount()= "+associationManager.getAssociatedAtoms().getAtomCount());
//    		throw new RuntimeException("rejectNotify wrong");
//    	}
    }

    public double energyChange() {return uNew - uOld;}
    
    public AtomIterator affectedAtoms() {
        return affectedAtomIterator;
    }

    public void setPressure(double p) {pressure = p;}
    public final double getPressure() {return pressure;}
    public Dimension getPressureDimension() {return Pressure.DIMENSION;}

	public Agent makeAgent(IAtom a, Box agentBox) {
		return new Agent();
	}

	public void releaseAgent(Agent agent, IAtom atom, Box agentBox) {
		
	}
	public static class Agent {
		public int nAtoms;//number of atoms in the smer
		public IAtom nextAtom;
	}
}
