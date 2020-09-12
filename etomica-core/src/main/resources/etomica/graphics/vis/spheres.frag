#version 330 core

layout (location = 0) out vec3 gPosition;
layout (location = 1) out vec3 gNormal;
layout (location = 2) out vec4 gAlbedoSpec;

struct Material {
//    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
//    float shininess;
};

uniform Material material;
uniform mat4 view;
uniform mat4 projection;
uniform vec3 viewPos;

flat in vec3 sphere_center;
flat in float radius_sq;
smooth in vec3 frag_pos;

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

    gPosition = frag_pos;
    gNormal = surface_normal;
    gAlbedoSpec.rgb = material.diffuse;
    gAlbedoSpec.a = material.specular.r;
}
