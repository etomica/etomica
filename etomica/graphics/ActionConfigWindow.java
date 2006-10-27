/**
 * 
 */
package etomica.graphics;

import java.awt.Color;
import java.awt.TextArea;

import javax.swing.JFrame;

import etomica.action.Action;
import etomica.atom.AtomLeaf;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.phase.Phase;
import etomica.space.Vector;

/**
 * Action that opens a new window and dumps the coordinates into the window.
 * @author andrew
 */
public class ActionConfigWindow implements Action {
    private final AtomIteratorLeafAtoms iterator;
    
    public ActionConfigWindow(Phase phase) {
        iterator = new AtomIteratorLeafAtoms(phase);
    }
    
    public String getLabel() {
        return "a label";
    }
    
    public void actionPerformed() {
        JFrame f = new JFrame();
        TextArea textArea = new TextArea();
        textArea.setEditable(false);
        textArea.setBackground(Color.white);
        textArea.setForeground(Color.black);
        iterator.reset();
        while (iterator.hasNext()) {
            AtomLeaf atom = (AtomLeaf)iterator.nextAtom();
            Vector pos = atom.coord.position();
            String str = Double.toString(pos.x(0));
            for (int i=1; i<pos.D(); i++) {
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