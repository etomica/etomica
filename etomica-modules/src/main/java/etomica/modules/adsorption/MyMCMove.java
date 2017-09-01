/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.adsorption;

import etomica.action.AtomActionRandomizeVelocity;
import etomica.box.Box;
import etomica.integrator.IntegratorBox;
import etomica.integrator.mcmove.MCMoveInsertDelete;
import etomica.math.DoubleRange;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeArrayList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

/**
 * @author kofke
 *
 * Extends insert/delete mcmove class by permitting insertion/deletion in a
 * subregion defined by a range of values of the z coordinate.
 */
public class MyMCMove extends MCMoveInsertDelete {

    private static final long serialVersionUID = 1L;
    private double zFraction, sigma;
    private Vector position;
    private final MoleculeArrayList activeAtoms;
    private final AtomActionRandomizeVelocity randomizer;
    private final IntegratorBox integrator;
    protected int testMoleculeIndex;
    protected DoubleRange range;
    protected final int dim;
    protected boolean bothSides;

	public MyMCMove(IntegratorBox integrator, IRandom random,
                    Space space, double zFraction, double sigma, int dim) {
		super(integrator.getPotentialMaster(), random, space);
		position = space.makeVector();
		setZFraction(zFraction, sigma);
		this.dim = dim;
        this.integrator = integrator;
        randomizer = new AtomActionRandomizeVelocity(0, random);
        activeAtoms = new MoleculeArrayList();
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
			    testMolecule = reservoir.remove(reservoir.getMoleculeCount()-1);
			}
			else {
			    testMolecule = species.makeMolecule();
			}
            position.E(positionSource.randomPosition());
			double z = position.getX(dim);
            double zBoundary = box.getBoundary().getBoxSize().getX(dim);
            z *= zFraction;
            if (bothSides) {
                if (z < 0) {
                    z += zBoundary*(0.5-zFraction) - sigma;
                }
                else {
                    z -= zBoundary*(0.5-zFraction) - sigma;
                }
            }
            else {
                z += 0.5*zBoundary*(1 - zFraction) - sigma;
            }
			position.setX(dim,z);
			atomTranslator.setDestination(position);
			atomTranslator.actionPerformed(testMolecule);
            box.addMolecule(testMolecule);
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

    public double getChi(double temperature) {//note that moleculeCount() gives the number of molecules after the trial is attempted
        if (insert) {
            energyMeter.setTarget(testMolecule);
            uNew = energyMeter.getDataAsScalar();
        }
        else {
            uNew = 0;
        }
        double b = uOld - uNew;
        if (insert) b += mu;
        else b -= mu;

        double a = insert ? zFraction * box.getBoundary().volume() / (activeAtoms.getMoleculeCount() + 1)
                : activeAtoms.getMoleculeCount() / zFraction / box.getBoundary().volume();
        return a * Math.exp(b / temperature);
    }

	public void rejectNotify() {
//        if(!insert) {
//            System.out.println("rejected deleting "+testMolecule.getIndex());
//        }
//        else {
//            System.out.println("rejected inserting "+testMolecule.getIndex());
//        }
	    super.rejectNotify();
	}
	
	public void acceptNotify() {
        super.acceptNotify();
		if(!insert) {
//	        System.out.println("accepted deleting "+testMolecule.getIndex());
			activeAtoms.remove(testMoleculeIndex);
		} else {
//            System.out.println("accepted inserting "+testMolecule.getIndex());
			activeAtoms.add(testMolecule);
            randomizer.setTemperature(integrator.getTemperature());
			randomizer.actionPerformed(testMolecule.getChildList().getAtom(0));
		}
	}
    
    public void setupActiveAtoms() {
    	activeAtoms.clear();
    	double zBoundary = box.getBoundary().getBoxSize().getX(dim);
        int nMolecules = moleculeList.getMoleculeCount();
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = moleculeList.getMolecule(i);
            if (molecule.getType() != species) continue;

    		double z = molecule.getChildList().getAtom(0).getPosition().getX(dim);
    		if (bothSides) {
    		    if (Math.abs(z) < 0.5*zBoundary*zFraction) continue;
    		}
    		else {
    		    if (z < zBoundary*(0.5-zFraction) - sigma || z > zBoundary*0.5 - sigma) continue;
    		}
    		activeAtoms.add(molecule);
    	}
    }

    /**
	 * Sets the zFraction, the fraction of the box volume into which atoms are
	 * inserted or deleted.  The volume is on the far left or right side (in
	 * the z dimension).  To put the volume on the left side (negative z),
	 * specify a negative z fraction (a positive value will result in the
	 * volume being on the positive z side.
	 */
	public void setZFraction(double zFraction, double zFractionMax) {
	    this.zFraction = zFraction;
	    this.sigma = zFractionMax;
	}
	
	public double getZFraction() {
	    return zFraction;
	}

	public void setSpecies(ISpecies s) {
		super.setSpecies(s);
		if (box != null) {
		    moleculeList = box.getMoleculeList(s);
		}
	}
}
