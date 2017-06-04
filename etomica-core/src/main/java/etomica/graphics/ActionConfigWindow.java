/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;
import java.awt.TextArea;

import javax.swing.JFrame;

import etomica.action.IAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Vector;
import etomica.box.Box;

/**
 * Action that opens a new window and dumps the coordinates into the window.
 * @author Andrew Schultz
 */
public class ActionConfigWindow implements IAction {
    private final IAtomList leafList;
    
    public ActionConfigWindow(Box box) {
        leafList = box.getLeafList();
    }
    
    public void actionPerformed() {
        JFrame f = new JFrame();
        TextArea textArea = new TextArea();
        textArea.setEditable(false);
        textArea.setBackground(Color.white);
        textArea.setForeground(Color.black);
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            Vector pos = a.getPosition();
            String str = Double.toString(pos.getX(0));
            for (int i=1; i<pos.getD(); i++) {
                str += " "+Double.toString(pos.getX(i));
            }
            textArea.append(str+"\n");
        }
        f.add(textArea);
        f.pack();
        f.setSize(400,600);
        f.setVisible(true);
    }
}
