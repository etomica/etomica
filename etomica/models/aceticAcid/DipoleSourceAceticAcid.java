/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.aceticAcid;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.DipoleSource;
import etomica.atom.MoleculeOriented;
import etomica.space.ISpace;
import etomica.space3d.IOrientationFull3D;
import etomica.units.Electron;

/**
 * dipole for acetic acid
 *
 * @author Hye Min Kim
 */
public class DipoleSourceAceticAcid implements DipoleSource {

    public DipoleSourceAceticAcid(ISpace space) {
        dipole = space.makeVector();
    }
    
    public IVector getDipole(IMolecule molecule) {//dipole= sum of position * charge for all sites in a molecule
        IAtomList childList = molecule.getChildList();

        IAtom cH3 = childList.getAtom(0);
        IAtom c = childList.getAtom(1);
        IAtom dBO = childList.getAtom(2);
        IAtom sBO = childList.getAtom(3);
        IAtom h = childList.getAtom(4);
		double zCH3 =  Electron.UNIT.toSim(0.08);//partial charge of CH3 site
		double zC =  Electron.UNIT.toSim(0.55);
        double zDBO   = Electron.UNIT.toSim(-0.50);
        double zSBO   = Electron.UNIT.toSim(-0.58);
        double zH   =  Electron.UNIT.toSim(0.45);
        
        dipole.Ea1Tv1(zCH3, cH3.getPosition());
        dipole.PEa1Tv1(zC, c.getPosition());
        dipole.PEa1Tv1(zDBO, dBO.getPosition());
        dipole.PEa1Tv1(zSBO, sBO.getPosition());
        dipole.PEa1Tv1(zH, h.getPosition());

        return dipole;
    }

    protected final IVectorMutable dipole;
}
