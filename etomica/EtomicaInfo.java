package etomica;

public class EtomicaInfo {
    
    private String description;
    
    public EtomicaInfo() {
        this("No description available");
    }
    public EtomicaInfo(String desc) {
        description = desc;
    }
    
    public String getDescription() {return description;}
    public void setDescription(String str) {description = str;}
    
}