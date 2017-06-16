/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.rotation;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;
import etomica.config.Configuration;
import etomica.space3d.Vector3D;

/**
 * reads configuration coordinates from a file and assigns them to the leaf atoms in a box
 */
public class ConfigurationWater3_3P implements Configuration, java.io.Serializable {

    public ConfigurationWater3_3P() {
    }
    
    public void initializeCoordinates(Box box) {
        Vector[] vec = new Vector[9];
        vec[ 0 ] = new Vector3D( -1.7018324526727353+.470824 , -2.2661777994405443 , -0.5630320704547286 );
        vec[ 1 ] = new Vector3D( -0.6908808660596318+.470824 , -1.0281664180026906 , -0.8990654429377792 );
        vec[ 2 ] = new Vector3D( -1.4882787859907576+.470824 , -1.3094875114788505 , -0.36519748205710423 );
        vec[ 3 ] = new Vector3D( -1.4876283712348002 , 0.34704039600242287 , 0.20855299296612284 );
        vec[ 4 ] = new Vector3D( -1.8292349266169308 , 1.8799977518289714 , 0.6568390902352385 );
        vec[ 5 ] = new Vector3D( -1.124213786523075 , 1.2754350591391652 , 0.2860977482885768 );
        vec[ 6 ] = new Vector3D( 0.14608749547957278 , 0.6039787851516225 , -0.6351451951992961 );
        vec[ 7 ] = new Vector3D( 1.3818907019750226 , 0.8142906064516727 , -1.6821280465892123 );
        vec[ 8 ] = new Vector3D( 0.7192645620836094 , 0.15747382037032204 , -1.3222415252323139 );
        
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            a.getPosition().E(vec[iLeaf]);
        }
    }
        
    private static final long serialVersionUID = 2L;
}
