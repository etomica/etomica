package simulate;
import simulate.beans.*;
import simulate.units.Unit;
import java.awt.*;
import java.beans.*;

public class DisplayBoxBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(Meter.class, MeterEditor.class);
        //PropertyEditorManager.registerEditor(Unit.class, DeviceUnitEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            return new PropertyDescriptor[] {
                new PropertyDescriptor("background", DisplayBox.class),
                new PropertyDescriptor("label",      DisplayBox.class),
                new PropertyDescriptor("meter",      DisplayBox.class),
                new PropertyDescriptor("name",       DisplayBox.class),
                new PropertyDescriptor("precision",  DisplayBox.class),
                new PropertyDescriptor("unit",       DisplayBox.class),
                new PropertyDescriptor("updateInterval", DisplayBox.class),
                new PropertyDescriptor("useCurrentValue", DisplayBox.class)};
        }
        catch(IntrospectionException e) {
            System.out.println("Error: "+e);
            return null;
        }
    }
    
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

