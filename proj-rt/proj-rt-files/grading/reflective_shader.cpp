#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& same_side_normal,int recursion_depth) const
{
    vec3 color;
    // TODO: determine the color
    Render_World render;
    Ray intersectRay;
    intersectRay.endpoint = intersection_point;
    intersectRay.direction = same_side_normal;
    color = render.Cast_Ray(intersectRay, recursion_depth) + reflectivity*(render.Cast_Ray(ray, recursion_depth) - render.Cast_Ray(intersectRay, recursion_depth));
    return color;
}
