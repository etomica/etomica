package etomica;
import etomica.*;
import java.awt.*;
import java.beans.*;

public class ConfigurationSequentialBeanInfo extends SimpleBeanInfo {
    
//    static {
//        PropertyEditorManager.registerEditor(ConfigurationSequential.class, ConfigurationSequentialEditor.class);
//    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
//            PropertyDescriptor boundaryDescriptor =
//                new PropertyDescriptor("boundary", Phase.class);
//            boundaryDescriptor.setPropertyEditorClass(BoundaryEditor.class);
            return new PropertyDescriptor[] {
                new PropertyDescriptor("fillVertical",ConfigurationSequential.class),
                };
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
    
//    public BeanDescriptor getBeanDescriptor() {
//        return new BeanDescriptor(Phase.class, PhaseCustomizer.class);
//    }
    
    
    public Image getIcon(int iconType){
        String name = "";
        if(iconType == BeanInfo.ICON_COLOR_16x16){
            name="COLOR_16x16";
        }
        else if (iconType == BeanInfo.ICON_COLOR_32x32){
            name="COLOR_32x32";
        }
        else if (iconType == BeanInfo.ICON_MONO_32x32){
            name="MONO_32x32";
        }
        else if (iconType == BeanInfo.ICON_MONO_16x16){
            name="MONO_16x16";
        }
        else return null;
        return loadImage("Molecule_" + name + ".gif");
    }
}

