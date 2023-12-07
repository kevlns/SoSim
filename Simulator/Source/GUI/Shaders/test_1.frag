#version 330 core
out vec4 FragColor;

void main()
{
    //    float dist = length(gl_PointCoord - vec2(0.5));
    //    if (dist > 0.5) discard;
    FragColor = vec4(0.3, 0.8, 0.8, 1.0);
}
