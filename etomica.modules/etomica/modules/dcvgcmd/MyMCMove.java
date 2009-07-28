package etomica.modules.dcvgcmd;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IRandom;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.atom.MoleculeArrayList;
import etomica.integrator.IntegratorBox;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.space.ISpace;

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
	public MyMCMove(IntegratorBox integrator, IRandom random,
			        ISpace space, double zFraction) {
		super(integrator.getPotentialMaster(), random, space);
		position = space.makeVector();
		setZFraction(zFraction);
        this.integrator = integrator;
        randomizer = new AtomActionRandomizeVelocity(0, random);
        activeAtoms = new MoleculeArrayList();
	}

    public void setBox(IBox p) {
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
			    testMolecule = reservoir.remove(reservoir.getMoleculeCount()-1);
			}
			else {
			    testMolecule = species.makeMolecule();
			}
            box.addMolecule(testMolecule);
			position.E(positionSource.randomPosition());
			double z = position.getX(2);
            z *= zFraction;
			if(nearOrigin) {
                z = (0.5*zFraction-0.5)*box.getBoundary().getBoxSize().getX(2) + z;
			} else {
				z = (0.5-0.5*zFraction)*box.getBoundary().getBoxSize().getX(2) - z;
			}
			position.setX(2,z); //multiply z-coordinate by zFraction
			atomTranslator.setDestination(position);
			atomTranslator.actionPerformed(testMolecule);
		} else {//delete
			if(activeAtoms.getMoleculeCount() == 0) {
				testMolecule = null;//added this line 09/19/02
				return false;
			}
            testMoleculeIndex = random.nextInt(activeAtoms.getMoleculeCount());
			testMolecule = activeAtoms.getMolecule(testMoleculeIndex);
			energyMeter.setTarget(testMolecule);
			uOld = energyMeter.getDataAsScalar();
		} 
		uNew = Double.NaN;
		return true;
	}//end of doTrial

	public double getA() {//note that moleculeCount() gives the number of molecules after the trial is attempted
		return insert ? zFraction*box.getBoundary().volume()/(activeAtoms.getMoleculeCount()+1) 
					  : activeAtoms.getMoleculeCount()/zFraction/box.getBoundary().volume();        
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
    	double zBoundary = box.getBoundary().getBoxSize().getX(2);
    	double zmin = nearOrigin ? -0.5*zBoundary : 0.5*(1.0-zFraction)*zBoundary;
    	double zmax = nearOrigin ? -0.5*(1.0-zFraction)*zBoundary : 0.5*zBoundary;
        int nMolecules = moleculeList.getMoleculeCount();
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = moleculeList.getMolecule(i);

    		double z = ((IAtomPositioned)molecule.getChildList().getAtom(0)).getPosition().getX(2);
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
	private IVectorMutable position;
	private boolean nearOrigin;
	private final MoleculeArrayList activeAtoms;
    private IMoleculeList moleculeList;
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
	public void setSpecies(ISpecies s) {
		super.setSpecies(s);
		moleculeList = box.getMoleculeList(s);
	}
}
