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
public class ConfigurationWater3 implements Configuration, java.io.Serializable {

    public ConfigurationWater3() {
    }
    
    public void initializeCoordinates(Box box) {
        Vector[] vec = new Vector[324];
        vec[ 0 ] = new Vector3D( -1.1533346461152898+.461317 , -0.400516942962945 , -0.31381683967094204 );
        vec[ 1 ] = new Vector3D( -1.3230564054265321+.461317 , -1.7697924565967533 , -0.936853182075963 );
        vec[ 2 ] = new Vector3D( -0.9934520861162452+.461317 , -0.8901236178710313 , -1.1206344838096907 );
        vec[ 3 ] = new Vector3D( -1.0561123122628668+.461317 , -0.940056280888558 , -0.9938258672915992 );
        vec[ 4 ] = new Vector3D( -1.5352592825444336 , 1.7594363991756543 , 1.0945909569886423 );
        vec[ 5 ] = new Vector3D( -0.2987710571614241 , 1.1800070379358236 , 0.44094451157596304 );
        vec[ 6 ] = new Vector3D( -1.0749175692641755 , 0.933795643034474 , 0.9441426641674294 );
        vec[ 7 ] = new Vector3D( -1.0344907448714007 , 1.0710056501236005 , 0.8989864262360243 );
        vec[ 8 ] = new Vector3D( 0.3628673472177766 , 0.42828824265580245 , -1.3502538711569134 );
        vec[ 9 ] = new Vector3D( 1.786816748275948 , 0.9336310596550622 , -1.2559568281938271 );
        vec[ 10 ] = new Vector3D( 0.8921778211692928 , 1.1357709735946988 , -0.9821092696910503 );
        vec[ 11 ] = new Vector3D( 0.9389442708130433 , 1.0193283005722875 , -1.064292011718281 );
        
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.size();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.get(iLeaf);
            a.getPosition().E(vec[iLeaf]);
        }
    }
        
    private static final long serialVersionUID = 2L;
}
