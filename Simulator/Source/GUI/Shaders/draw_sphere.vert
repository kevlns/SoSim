#version 330 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec4 aColor;

out vec3 pointCenter;
out vec4 pointColor;

uniform vec4 pureColor;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float particleRadius;

void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0f);
    pointCenter = vec3(gl_Position);
    pointColor = pureColor;
}