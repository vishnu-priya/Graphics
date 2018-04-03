#include "plane.h"
#include "ray.h"
#include <cfloat>


// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
bool Plane::
Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // TODO
    vec3 del = x1 - ray.endpoint;
    Hit hit;
    hit.object = this;
    hit.ray_exiting = 0;
    double t = dot(del, normal);
    double n = dot(ray.direction, normal);
    if(n){
      if(t/n > small_t){
        hit.t = t/n;
      }else{
        return false;
      }
      hits.push_back(hit);
      return true;
    }
    return false;
}

vec3 Plane::
Normal(const vec3& point) const
{
    return normal;
}
