package simulate.gui;

import simulate.*;
import java.awt.Panel;

public class FormDesignFrame extends javax.swing.JInternalFrame {
    
    public FormDesignFrame(String space){
        setResizable(true);
        setIconifiable(true);
        setMaximizable(true);
        setClosable(true);
        setBounds(230, 200, 775, 400);
        setTitle("Form Design");
    }
}