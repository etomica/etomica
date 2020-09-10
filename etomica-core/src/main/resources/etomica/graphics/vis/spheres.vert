#version 330 core

layout (location = 0) in vec3 position;

uniform mat4 view;

out float vRadius;

void main() {
    gl_Position = view * vec4(position, 1);
    vRadius = 1;
}
