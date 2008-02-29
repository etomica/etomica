package etomica.graphics;

import java.awt.Color;
import java.awt.TextArea;

import javax.swing.JFrame;

import etomica.action.Action;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;

/**
 * Action that opens a new window and dumps the coordinates into the window.
 * @author Andrew Schultz
 */
public class ActionConfigWindow implements Action {
    private final AtomSet leafList;
    
    public ActionConfigWindow(IBox box) {
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
            IAtomPositioned a = (IAtomPositioned)leafList.getAtom(iLeaf);
            IVector pos = a.getPosition();
            String str = Double.toString(pos.x(0));
            for (int i=1; i<pos.getD(); i++) {
                str += " "+Double.toString(pos.x(i));
            }
            textArea.append(str+"\n");
        }
        f.add(textArea);
        f.pack();
        f.setSize(400,600);
        f.setVisible(true);
    }
}