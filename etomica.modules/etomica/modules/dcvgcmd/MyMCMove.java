package etomica.modules.dcvgcmd;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomSet;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.util.IRandom;

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
	public MyMCMove(IntegratorBox integrator, IRandom random, double zFraction) {
		super(integrator.getPotential(), random);
		position =  (Vector3D)integrator.getPotential().getSpace().makeVector();
		setZFraction(zFraction);
        this.integrator = integrator;
        randomizer = new AtomActionRandomizeVelocity(0, random);
        activeAtoms = new AtomArrayList();
	}

    public void setBox(Box p) {
        super.setBox(p);
        energyMeter.setBox(p);
    }
    
	/**
	 * Chooses and performs with equal probability an elementary molecule insertion
	 * or deletion.
	 */
	public boolean doTrial() {
		insert = (random.nextInt(2) == 0);
		if(insert) {
			uOld = 0.0;
			if(!reservoir.isEmpty()) {
			    testMolecule = (IMolecule)reservoir.remove(reservoir.getAtomCount()-1);
			}
			else {
			    testMolecule = species.makeMolecule();
			}
            box.addMolecule(testMolecule);
			position = (Vector3D)box.getBoundary().randomPosition();
			double z = position.x(2);
            z *= zFraction;
			if(nearOrigin) {
                z = (0.5*zFraction-0.5)*box.getBoundary().getDimensions().x(2) + z;
			} else {
				z = (0.5-0.5*zFraction)*box.getBoundary().getDimensions().x(2) - z;
			}
			position.setX(2,z); //multiply z-coordinate by zFraction
			atomTranslator.setDestination(position);
			atomTranslator.actionPerformed(testMolecule);
		} else {//delete
			if(activeAtoms.getAtomCount() == 0) {
				testMolecule = null;//added this line 09/19/02
				return false;
			}
            testMoleculeIndex = random.nextInt(activeAtoms.getAtomCount());
			testMolecule = (IMolecule)activeAtoms.getAtom(testMoleculeIndex);
			energyMeter.setTarget(testMolecule);
			uOld = energyMeter.getDataAsScalar();
		} 
		uNew = Double.NaN;
		return true;
	}//end of doTrial

	public double getA() {//note that moleculeCount() gives the number of molecules after the trial is attempted
		return insert ? zFraction*box.volume()/(activeAtoms.getAtomCount()+1) 
					  : activeAtoms.getAtomCount()/zFraction/box.volume();        
	}

	public void acceptNotify() {
        super.acceptNotify();
		if(!insert) {
			activeAtoms.remove(testMoleculeIndex);
			deltaN--;
		} else {
			activeAtoms.add(testMolecule);
            randomizer.setTemperature(integrator.getTemperature());
			randomizer.actionPerformed(testMolecule.getChildList().getAtom(0));
			deltaN++;
		}
	}
    
    public void setupActiveAtoms() {
    	activeAtoms.clear();
    	double zBoundary = box.getBoundary().getDimensions().x(2);
    	double zmin = nearOrigin ? -0.5*zBoundary : 0.5*(1.0-zFraction)*zBoundary;
    	double zmax = nearOrigin ? -0.5*(1.0-zFraction)*zBoundary : 0.5*zBoundary;
        int nMolecules = moleculeList.getAtomCount();
        for (int i=0; i<nMolecules; i++) {
            IAtom molecule = moleculeList.getAtom(i);

    		double z = ((IAtomPositioned)((IMolecule)molecule).getChildList().getAtom(0)).getPosition().x(2);
    		if(z < zmin || z > zmax) continue;
    		activeAtoms.add(molecule);
    	}
    }
    
    public int getDeltaN() {
    	return deltaN;
    }

    private static final long serialVersionUID = 1L;
	private double zFraction;
	private int deltaN = 0;
	private Vector3D position;
	private boolean nearOrigin;
	private final AtomArrayList activeAtoms;
    private AtomSet moleculeList;
	private final AtomActionRandomizeVelocity randomizer;
    private final IntegratorBox integrator;
    protected int testMoleculeIndex;
	
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
		moleculeList = box.getMoleculeList(s);
	}
}
