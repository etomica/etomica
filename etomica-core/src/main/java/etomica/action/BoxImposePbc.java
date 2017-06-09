/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.space.Vector;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.atom.IMoleculePositionDefinition;
import etomica.atom.IMoleculePositioned;
import etomica.space.Space;

/**
 * Action that imposes the central-image effect of a box having periodic
 * boundaries. Causes all atoms with coordinates outside the box boundaries to
 * be moved to the central-image location (inside the boundaries).
 */

public class BoxImposePbc extends BoxActionAdapter {

	/**
	 * Creates the action without specifying a box. Requires call to setBox
	 * before action can have any effect. Default is to apply central-imaging at
	 * the atom rather than molecule level.
	 */
	public BoxImposePbc(Space space) {
		setApplyToMolecules(false);
		this.space = space;
		setPositionDefinition(new MoleculePositionGeometricCenter(space));
	}
    
    public int getPriority() {return 100;}//100-199 is priority range for classes imposing PBC

	/**
	 * Creates the action ready to perform on the given box.
	 * 
	 * @param box
	 */
	public BoxImposePbc(Box box, Space space) {
		this(space);
		setBox(box);
	}

	public void actionPerformed() {
		Boundary boundary = box.getBoundary();
        if (applyToMolecules) {
            IMoleculeList molecules = box.getMoleculeList();
            
            for (int i=0; i<molecules.getMoleculeCount(); i++) {
                IMolecule molecule = molecules.getMolecule(i);
                Vector shift;
                if (molecule instanceof IMoleculePositioned) {
                    Vector position = ((IMoleculePositioned)molecule).getPosition();
                    shift = boundary.centralImage(position);
                    position.PE(shift);
                }
                else {
                    shift = boundary.centralImage(positionDefinition.position(molecule));
                }
                if (!shift.isZero()) {
                    translator.setTranslationVector(shift);
                    moleculeTranslator.actionPerformed(molecule);
                }
            }
        }
        else {
            IAtomList atoms = box.getLeafList();
            for (int i=0; i<atoms.getAtomCount(); i++) {
                IAtom atom = atoms.getAtom(i);
                Vector shift = boundary.centralImage(atom.getPosition());
                if (!shift.isZero()) {
                    atom.getPosition().PE(shift);
                }
            }
        }
	}

	/**
	 * Returns the value of applyToMolecules.
	 * 
	 * @return boolean
	 */
	public boolean isApplyToMolecules() {
		return applyToMolecules;
	}

	/**
	 * Sets a flag indicating whether periodic boundaries are applied to the
	 * molecules (true), or to the atoms (false). If applied to the atoms (the
	 * default case), then central imaging is done to each atom individually,
	 * which could cause a molecule to be split, with some of its atoms on one
	 * edge of the simulation box, and others on the other edge. If applied to
	 * molecules, the entire molecule will be shifted as a whole when enforcing
	 * central imaging.
	 * 
	 * @param applyToMolecules
	 *            The new value of the flag.
	 */
	public void setApplyToMolecules(boolean applyToMolecules) {
		this.applyToMolecules = applyToMolecules;
		if (applyToMolecules) {
	        translator = new AtomActionTranslateBy(space);
	        moleculeTranslator = new MoleculeChildAtomAction(translator);
		}
	}

    public void setPositionDefinition(IMoleculePositionDefinition positionDefinition) {
        this.positionDefinition = positionDefinition;
    }

    public IMoleculePositionDefinition getPositionDefinition() {
        return positionDefinition;
    }

    private static final long serialVersionUID = 1L;
    private AtomActionTranslateBy translator;
    private MoleculeChildAtomAction moleculeTranslator;
    private Space space;
    private IMoleculePositionDefinition positionDefinition;

	private boolean applyToMolecules;
}
