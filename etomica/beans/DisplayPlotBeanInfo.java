package simulate;
import simulate.beans.*;
import simulate.units.Unit;
import java.awt.*;
import java.beans.*;

public class DisplayPlotBeanInfo extends SimpleBeanInfo {
    
    static {
        PropertyEditorManager.registerEditor(MeterFunction.class, MeterFunctionEditor.class);
    }
    
    public PropertyDescriptor[] getPropertyDescriptors() {
        try {
            return new PropertyDescriptor[] {
                new PropertyDescriptor("name",       DisplayPlot.class),
                new PropertyDescriptor("label",  DisplayPlot.class),
                new PropertyDescriptor("meter", DisplayPlot.class),
                new PropertyDescriptor("makeResetButton",      DisplayPlot.class),
                new PropertyDescriptor("useCurrentValue",      DisplayPlot.class),
                new PropertyDescriptor("connected", DisplayPlot.class),
                new PropertyDescriptor("updateInterval", DisplayPlot.class)};
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

