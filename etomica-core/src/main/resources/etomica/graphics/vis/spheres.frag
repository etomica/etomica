#version 330 core

out vec4 FragColor;


struct DirLight {
    vec3 direction;

    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};


struct Material {
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
    float shininess;
};

uniform DirLight dirLight;
uniform Material material;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 viewPos;

flat in vec3 sphere_center;
flat in float radius_sq;
smooth in vec3 frag_pos;

vec3 calcDirLight(DirLight light, vec3 normal, vec3 viewDir) {
    vec3 lightDir = normalize(-light.direction);
    lightDir = (view * vec4(lightDir, 0)).xyz;
    vec3 halfwayDir = normalize(lightDir + viewDir);

    float diffuseFac = max(dot(normal, lightDir), 0.0);

    vec3 reflectDir = reflect(-lightDir, normal);
    float specularFac = pow(max(dot(normal, halfwayDir), 0.0), material.shininess);

    vec3 ambient = light.ambient * material.diffuse;
    vec3 diffuse = light.diffuse * diffuseFac * material.diffuse;
    vec3 specular = light.specular * specularFac * material.specular;
    return (ambient + diffuse + specular);
}

void main() {
    vec3 sphere_dir = sphere_center;
    vec3 ray_dir = normalize(frag_pos);
    gl_FragDepth = gl_FragCoord.z;

    float b = dot(ray_dir, sphere_dir);
    float disc = b*b + radius_sq - dot(sphere_dir, sphere_dir);
    if (disc < 0.0) {
//        FragColor = vec4(1);
//        return;
        discard;
    }

    float t = b - sqrt(disc);
    vec3 hit = t * ray_dir;

    vec4 pos = vec4(hit, 1.0);
    vec4 screen_pos = projection * pos;
    gl_FragDepth = (screen_pos.z / screen_pos.w + 1.0) / 2.0;

    vec3 surface_normal = normalize(hit - sphere_center);

    vec3 viewDir = normalize(-frag_pos);
    vec3 color = calcDirLight(dirLight, surface_normal, viewDir);
    float edgeFac = dot(viewDir, surface_normal);
    float fac = smoothstep(0.3, 0.4, edgeFac);
    color = color * fac;
    FragColor = vec4(color, 1.0);
}
