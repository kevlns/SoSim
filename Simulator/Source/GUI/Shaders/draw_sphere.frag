#version 330 core
out vec4 FragColor;

in vec3 pointCenter;
in vec4 pointColor;

uniform float particleRadius;

out vec4 color;

void main() {
    vec3 FragRayOrigin = vec3(5, 5, 5);
    vec3 FragRayDirection = vec3(-1, -1, -1);

    // 计算光线和球体的交点
    vec3 oc = FragRayOrigin - pointCenter;
    float a = dot(FragRayDirection, FragRayDirection);
    float b = 2.0 * dot(oc, FragRayDirection);
    float c = dot(oc, oc) - particleRadius * particleRadius;
    float discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {
        // 计算最近的交点（t值）
        float t = (-b - sqrt(discriminant)) / (2.0 * a);
        vec3 hitPoint = FragRayOrigin + t * FragRayDirection;

        // 计算法线，简单的光照和颜色
        vec3 normal = normalize(hitPoint - pointCenter);
        float lightIntensity = dot(normal, normalize(vec3(1.0, 1.0, 1.0))); // 光源方向
        color = pointColor + vec4(0.5 * lightIntensity + 0.5); // 简单的漫反射
    } else {
        // 如果没有交点，则显示背景颜色
        color = vec4(0.0, 0.0, 0.0, 1.0);
    }
}