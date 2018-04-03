#include "sphere.h"
#include "ray.h"
#include<math.h>


// Determine if the ray intersects with the sphere
bool Sphere::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    // TODO
    vec3 y = ray.endpoint - center;
    Hit h1, h2;
    h1.object = this;
    h1.ray_exiting = 0;
    h2.object = this;
    h2.ray_exiting = 1;
    double dott = dot(ray.direction, y);
    double delta = pow(dott, 2.0) - (dot(y, y) - pow(radius, 2.0));
    if(!(delta<0)){
      double t0 = -dott - delta;
      double t1 = -dott + delta;
      if(t0 > 0 || t1 > 0){
        if(t0 > 0 && t1 > 0){
          if(t0 > t1){
            h1.t = t1;
            h2.t = t0;
          }else{
            h1.t = t0;
            h2.t = t1;
          }
        }else if(t0 > 0){
          h1.t = 0;
          h2.t = t0;
        }else{
          h1.t = 0;
          h2.t = t1;
        }
        hits.push_back(h1);
        hits.push_back(h2);
        return true;
      }
    }
    return false;
}

vec3 Sphere::Normal(const vec3& point) const
{
    vec3 normal;
    // TODO: set the normal
    normal = cross(point, center);
    return normal;
}
