package etomica.modules.dcvgcmd;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.integrator.IntegratorPhase;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space3d.Vector3D;
import etomica.species.Species;

/**
 * @author kofke
 *
 * Extends insert/delete mcmove class by permitting insertion/deletion in a
 * subregion defined by a range of values of the z coordinate.
 */
public class MyMCMove extends MCMoveInsertDelete {

	/**
	 * Constructor for MyMCMove.
	 * @param parent
	 */
	public MyMCMove(IntegratorPhase integrator, double zFraction) {
		super(integrator.getPotential());
		position =  (Vector3D)integrator.getPotential().getSpace().makeVector();
		setZFraction(zFraction);
        this.integrator = integrator;
	}

    public void setPhase(Phase p) {
        super.setPhase(p);
        energyMeter.setPhase(p);
    }
    
	/**
	 * Chooses and performs with equal probability an elementary molecule insertion
	 * or deletion.
	 */
	public boolean doTrial() {
		insert = Simulation.random.nextDouble() < 0.5;
		if(insert) {
			uOld = 0.0;
			if(!reservoir.isEmpty()) testMolecule = reservoir.removeFirst();
			else testMolecule = moleculeFactory.makeAtom();
            phase.addMolecule(testMolecule, speciesAgent);
			position = (Vector3D)phase.randomPosition();
			double z = position.x(2);
			if(nearOrigin) {
				z *= zFraction;
			} else {
				z *= zFraction;
				z = phase.getBoundary().getDimensions().x(2) - z;
			}
			position.setX(2,z); //multiply z-coordinate by zFraction		
			atomTranslator.setDestination(position);
			atomTranslator.actionPerformed(testMolecule);
            phase.getCellManager().assignCell(testMolecule);
		} else {//delete
			if(activeAtoms.size() == 0) {
				testMolecule = null;//added this line 09/19/02
				return false;
			}
			testMolecule = activeAtoms.getRandom();
			energyMeter.setTarget(testMolecule);
			uOld = energyMeter.getDataAsScalar();
		} 
		uNew = Double.NaN;
		return true;
	}//end of doTrial

	public double getA() {//note that moleculeCount() gives the number of molecules after the trial is attempted
		return insert ? zFraction*phase.volume()/(activeAtoms.size()+1) 
					  : activeAtoms.size()/zFraction/phase.volume();        
	}

	public void acceptNotify() {
		if(!insert) {
			phase.removeMolecule(testMolecule);
			activeAtoms.remove(testMolecule);
			deltaN--;
		} else {
			activeAtoms.add(testMolecule);
            randomizer.setTemperature(integrator.getTemperature());
			randomizer.actionPerformed(testMolecule);
			deltaN++;
		}
	}
    
    public void setupActiveAtoms() {
    	activeAtoms.clear();
    	atomIterator.reset();
    	double zBoundary = phase.getBoundary().getDimensions().x(2);
    	double zmin = nearOrigin ? 0.0 : (1.0-zFraction)*zBoundary;
    	double zmax = nearOrigin ? zFraction*zBoundary : zBoundary;
    	while(atomIterator.hasNext()) {
    		Atom atom = atomIterator.nextAtom();
    		double z = atom.coord.position().x(2);
    		if(z < zmin || z > zmax) continue;
    		activeAtoms.add(atom);
    	}
    }
    
    public int getDeltaN() {
    	return deltaN;
    }

	private double zFraction;
	private int deltaN = 0;
	private Vector3D position;
	private boolean nearOrigin;
	private AtomList activeAtoms = new AtomList();
	private AtomIteratorArrayListSimple atomIterator;// = new AtomIteratorList();
	private final AtomActionRandomizeVelocity randomizer = new AtomActionRandomizeVelocity(0);
    private final IntegratorPhase integrator;
	
	/**
	 * Returns the zFraction.
	 * @return double
	 */
	public double getZFraction() {
		return nearOrigin ? -zFraction : +zFraction;
	}

	/**
	 * Sets the zFraction.
	 * @param zFraction The zFraction to set
	 */
	public void setZFraction(double zFraction) {
		this.zFraction = zFraction;
		nearOrigin = zFraction < 0.0;
		this.zFraction = Math.abs(zFraction);		
	}

	/**
	 * @see etomica.MCMoveInsertDelete#setSpecies(etomica.Species)
	 */
	public void setSpecies(Species s) {
		super.setSpecies(s);
		atomIterator = new AtomIteratorArrayListSimple();
		atomIterator.setList(((AtomTreeNodeGroup)speciesAgent.node).childList);
	}

}
