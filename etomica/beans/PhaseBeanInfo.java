package etomica;
import etomica.beans.*;
import java.awt.*;
import java.beans.*;

public class PhaseBeanInfo extends SimpleBeanInfo implements java.io.Serializable {
    
    static {
        PropertyEditorManager.registerEditor(Space.Boundary.class, BoundaryEditor.class);
        PropertyEditorManager.registerEditor(Space.Coordinate[].class, MoleculePositionEditor.class);
//        PropertyEditorManager.registerEditor(Configuration.class, ConfigurationEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            return new PropertyDescriptor[] {
                new PropertyDescriptor("name",       Phase.class),
                new PropertyDescriptor("boundary", Phase.class),
                new PropertyDescriptor("configuration", Phase.class),
                new PropertyDescriptor("density",      Phase.class)};
        }
/*        try {
            PropertyDescriptor boundaryDescriptor =
                new PropertyDescriptor("boundary", Phase.class);
            boundaryDescriptor.setPropertyEditorClass(BoundaryEditor.class);
//            IndexedPropertyDescriptor moleculePositionDescriptor =
//                new IndexedPropertyDescriptor("moleculePosition", Phase.class);
//            moleculePositionDescriptor.setPropertyEditorClass(MoleculePositionEditor.class);
//            PropertyDescriptor meterListDescriptor =
//                new PropertyDescriptor("meterList", Phase.class);
//            meterListDescriptor.setPropertyEditorClass(MeterListEditor.class);
//            return new PropertyDescriptor[] {boundaryDescriptor, meterListDescriptor};
            return new PropertyDescriptor[] {boundaryDescriptor,
//                                             new IndexedPropertyDescriptor("moleculePosition", Phase.class),
                                             //moleculePositionDescriptor,
                                             new PropertyDescriptor("registered", Phase.class)};
        }*/
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
    
    public BeanDescriptor getBeanDescriptor() {
        return new BeanDescriptor(Phase.class, PhaseCustomizer.class);
    }
    
}

