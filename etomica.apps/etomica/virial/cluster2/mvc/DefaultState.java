package etomica.virial.cluster2.mvc;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class DefaultState implements State {

  private Map<String, Object> data = new HashMap<String, Object>();

  private Map<String, Object> getData() {

    return data;
  }

  public void clear() {

    getData().clear();
  }

  public Object getProperty(String key) {

    return getData().get(key);
  }

  public Set<String> getKeys() {

    return getData().keySet();
  }

  public void setProperty(String key, Object value) {

    getData().put(key, value);
  }
}
