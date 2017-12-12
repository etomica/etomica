package etomica.server.representations;

import com.fasterxml.jackson.annotation.JsonProperty;

public class PropertyUpdate {
    public @JsonProperty long id;
    public @JsonProperty String property;
    public @JsonProperty Object newValue;


    public PropertyUpdate() {
        super();
    }

    public PropertyUpdate(long id, String property, Object newValue) {
        this.id = id;
        this.property = property;
        this.newValue = newValue;
    }

    public long getId() {
        return id;
    }

    public String getProperty() {
        return property;
    }

    public Object getNewValue() {
        return newValue;
    }
}
