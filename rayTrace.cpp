#include <cstdlib> 
#include <cstdio> 
#include <cmath> 
#include <fstream> 
#include <vector> 
#include <iostream> 
#include <cassert> 
 
#if defined __linux__ || defined __APPLE__ 
// "Compiled for Linux
#else 
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793 
#define INFINITY 1e8 
#endif 
 
template<typename T> 
class Vec3 
{ 
public: 
    T x, y, z; 
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {} 
    Vec3(T xx) : x(xx), y(xx), z(xx) {} 
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {} 
    Vec3& normalize() 
    { 
        T nor2 = length2(); 
        if (nor2 > 0) { 
            T invNor = 1 / sqrt(nor2); 
            x *= invNor, y *= invNor, z *= invNor; 
        } 
        return *this; 
    } 
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); } 
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); } 
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); } 
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); } 
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; } 
    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; } 
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); } 
    T length2() const { return x * x + y * y + z * z; } 
    T length() const { return sqrt(length2()); } 
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v) 
    { 
        os << "[" << v.x << " " << v.y << " " << v.z << "]"; 
        return os; 
    } 
}; 
 
typedef Vec3<float> Vec3f; 

/*
class Shape
{
    public:
    virtual bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const = 0;
    virtual Vec3f normal(const Vec3f &phit) const;
    virtual Vec3f lightDirection(const Vec3f &phit) const;
};*/


class Sphere
{ 
public: 
    Vec3f center;                           /// position of the sphere 
    float radius, radius2;                  /// sphere radius and radius^2 
    // Not sure how these will mesh with base class, could cause error
    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light) 
    float transparency, reflection;         /// surface transparency and reflectivity 
    Sphere( 
        const Vec3f &c, 
        const float &r, 
        const Vec3f &sc, 
        const float &refl = 0, 
        const float &transp = 0, 
        const Vec3f &ec = 0) : 
        center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec), 
        transparency(transp), reflection(refl) 
    { /* empty */ } 
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const 
    { 
        Vec3f l = center - rayorig; 
        float tca = l.dot(raydir); 
        if (tca < 0) return false; 
        float d2 = l.dot(l) - tca * tca; 
        if (d2 > radius2) return false; 
        float thc = sqrt(radius2 - d2); 
        t0 = tca - thc; 
        t1 = tca + thc; 
 
        return true; 
    }
    Vec3f normal(const Vec3f &phit) const{
        return phit - center;
    }
    Vec3f lightDirection(const Vec3f &phit) const{
        return center - phit;
    }
};


class Box{
    public:
    // define a center and 3 vectors
    // work from there
    Vec3f bounds [2];
    Vec3f surfaceColor;
    float transparency, reflection;
    Box(
        const Vec3f &b0, const Vec3f &b1,
        const Vec3f &sc, const float transp, const float refl) :
        surfaceColor(sc), transparency(transp), reflection(refl)
        { bounds[0] = b0, bounds[1] = b1; }
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t) const
    {
        float tmin, tmax, tymin, tymax, tzmin, tzmax;

        int sign0, sign1, sign2;

        Vec3f invdir = Vec3f(1/raydir.x, 1/raydir.y, 1/raydir.z);

        sign0 = invdir.x < 0;
        sign1 = invdir.y < 0;
        sign2 = invdir.z < 0;

        tmin = (bounds[sign0].x - rayorig.x) * invdir.x;
        tmax = (bounds[1-sign0].x -rayorig.x) * invdir.x;
        tymin = (bounds[sign1].y - rayorig.y) * invdir.y;
        tymax = (bounds[1-sign1].y - rayorig.y) *invdir.y;

        if ((tmin > tymax) || (tymin > tmax)){
            return false;
        }

        if(tymin > tmin) tmin = tymin;
        if(tymax < tmax) tmax = tymax;

        tzmin = (bounds[sign2].z - rayorig.z) * invdir.z;
        tzmax = (bounds[1-sign2].z - rayorig.z) * invdir.z;

        if((tmin > tzmax) || (tzmin > tmax)) return false;
        if(tzmin > tmin) tmin = tzmin;
        if(tzmax < tmax) tmax = tzmax;

        t = tmin;

        if(t < 0) {
            t = tmax;
            if (t < 0) return false;
        }

        return true;
    }
    Vec3f normal(const Vec3f &phit) const{
        // take the point hit, figure out what face its on, then return the normal of that face
        /*
        Ok, assume that we have two points (0,0,0) and (1,1,1) which define a box in 3d space
        We need to know what face of the box we have hit. the 8 points of our box can be defined
        as (0,0,1)(0,1,0)(1,0,0)(1,1,0)(1,0,1)(0,1,1). We take each point and create three others by
        replacing one of its values with one of its counterparts. We can use these points to calculate
        the normals on each plane. Then with phit, we can determine if it is on the plane, by checking if
        the dot product of (phit - point on plane) and the normal is 0, if it is return that normal
        */
        // we are going to check as we calculate each plane and normal

        // This is the naive implementation but I don't have time for more
        // Also not sure if my value for epsilon is sufficient there could be
        // edge cases (pobably litterally) where weird stuff happens

        /*
        So we have a vector which is parrallel to the normal that we want. Vectors are parallel when they
        are scalar multiples of each other, so our question is norm*c = phit. Therefore phit/norm = c
        */

        //1
        Vec3f pointOnPlane = Vec3f(bounds[0].x, bounds[0].y, bounds[1].z);
        Vec3f norm = pointOnPlane - bounds[0];
        float d = (phit - pointOnPlane).dot(norm);
        if(abs(d) < 0.0000001) {
            float scalar = phit.x/norm.x;
            return norm * scalar;
        }

        //2
        pointOnPlane = Vec3f(bounds[0].x, bounds[1].y, bounds[0].z);
        norm = pointOnPlane - bounds[0];
        d = (phit - pointOnPlane).dot(norm);
        if(abs(d) < 0.0000001) {
            float scalar = phit.x/norm.x;
            return norm * scalar;
        }
        //3
        pointOnPlane = Vec3f(bounds[1].x, bounds[0].y, bounds[0].z);
        norm = pointOnPlane - bounds[0];
        d = (phit - pointOnPlane).dot(norm);
        if(abs(d) < 0.0000001) {
            float scalar = phit.x/norm.x;
            return norm * scalar;
        }

        //4
        pointOnPlane = Vec3f(bounds[1].x, bounds[1].y, bounds[0].z);
        norm = pointOnPlane - bounds[1];
        d = (phit - pointOnPlane).dot(norm);
        if(abs(d) < 0.0000001) {
            float scalar = phit.x/norm.x;
            return norm * scalar;
        }

        //5
        pointOnPlane = Vec3f(bounds[1].x, bounds[0].y, bounds[1].z);
        norm = pointOnPlane - bounds[1];
        d = (phit - pointOnPlane).dot(norm);
        if(abs(d) < 0.0000001) {
            float scalar = phit.x/norm.x;
            return norm * scalar;
        }

        //6
        pointOnPlane = Vec3f(bounds[0].x, bounds[1].y, bounds[1].z);
        norm = pointOnPlane - bounds[1];
        d = (phit - pointOnPlane).dot(norm);
        if(abs(d) < 0.0000001) {
            float scalar = phit.x/norm.x;
            return norm * scalar;
        }

        std::cerr << "Attemt to calculate normal on box failed" << std::endl;
        throw std::exception();
    }
    // We will not allow boxes to be lights
};

#define MAX_RAY_DEPTH 5 
 
float mix(const float &a, const float &b, const float &mix) 
{ 
    return b * mix + a * (1 - mix); 
} 

Vec3f trace( 
    const Vec3f &rayorig, 
    const Vec3f &raydir, 
    const std::vector<Sphere> &spheres,
    const std::vector<Box> &boxes, 
    const int &depth) 
{ 
    //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    float tnear = INFINITY; 
    const Sphere* sphere = NULL; 
    const Box* box = NULL;
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < spheres.size(); ++i) { 
        float t0 = INFINITY, t1 = INFINITY; 
        if (spheres[i].intersect(rayorig, raydir, t0, t1)) { 
            if (t0 < 0) t0 = t1; 
            if (t0 < tnear) { 
                tnear = t0; 
                sphere = &spheres[i]; 
            } 
        }
    }

    for(unsigned i = 0; i < boxes.size(); ++i) {
        float t = INFINITY;
        if(boxes[i].intersect(rayorig, raydir, t)){
            if(t < tnear){
                tnear = t;
                box = &boxes[i];
            }
        }
    }
    // if there's no intersection return black or background color
    if (!sphere && !box) return Vec3f(2);
    else if (!box){
        Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray 
        Vec3f phit = rayorig + raydir * tnear; // point of intersection 
        Vec3f nhit = sphere->normal(phit); // normal at the intersection point 
        nhit.normalize(); // normalize normal direction 
        // If the normal and the view direction are not opposite to each other
        // reverse the normal direction. That also means we are inside the sphere so set
        // the inside bool to true. Finally reverse the sign of IdotN which we want
        // positive.
        float bias = 1e-4; // add some bias to the point from which we will be tracing 
        bool inside = false; 
        if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true; 
        if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) { 
            float facingratio = -raydir.dot(nhit); 
            // change the mix value to tweak the effect
            float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1); 
            // compute reflection direction (not need to normalize because all vectors
            // are already normalized)
            Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit); 
            refldir.normalize(); 
            Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, boxes, depth + 1); 
            Vec3f refraction = 0; 
            // if the sphere is also transparent compute refraction ray (transmission)
            if (sphere->transparency) { 
                float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface? 
                float cosi = -nhit.dot(raydir); 
                float k = 1 - eta * eta * (1 - cosi * cosi); 
                Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k)); 
                refrdir.normalize(); 
                refraction = trace(phit - nhit * bias, refrdir, spheres, boxes, depth + 1); 
            } 
            // the result is a mix of reflection and refraction (if the sphere is transparent)
            surfaceColor = ( 
                reflection * fresneleffect + 
                refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor; 
        } 
        else { 
            // it's a diffuse object, no need to raytrace any further
            for (unsigned i = 0; i < spheres.size(); ++i) { 
                if (spheres[i].emissionColor.x > 0) { 
                    // this is a light
                    Vec3f transmission = 1; 
                    Vec3f lightDirection = spheres[i].lightDirection(phit); 
                    lightDirection.normalize(); 
                    for (unsigned j = 0; j < spheres.size(); ++j) { 
                        if (i != j) { 
                            float t0, t1; 
                            if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) { 
                                transmission = 0; 
                                break; 
                            } 
                        } 
                    } 
                    surfaceColor += sphere->surfaceColor * transmission * 
                    std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor; 
                } 
            } 
        } 
    
        return surfaceColor + sphere->emissionColor; 
    } else {
        // in this case you have hit a box not a sphere
        Vec3f surfaceColor = 0;
        Vec3f phit = rayorig + raydir * tnear;
        Vec3f nhit = box->normal(phit);
        nhit.normalize();

        float bias = 1e-4;
        bool inside = false;
        if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true; 
        if ((box->transparency > 0 || box->reflection > 0) && depth < MAX_RAY_DEPTH) { 
            float facingratio = -raydir.dot(nhit); 
            // change the mix value to tweak the effect
            float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1); 
            // compute reflection direction (not need to normalize because all vectors
            // are already normalized)
            Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit); 
            refldir.normalize(); 
            Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, boxes, depth + 1); 
            Vec3f refraction = 0; 
            // if the sphere is also transparent compute refraction ray (transmission)
            if (box->transparency) { 
                float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface? 
                float cosi = -nhit.dot(raydir); 
                float k = 1 - eta * eta * (1 - cosi * cosi); 
                Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k)); 
                refrdir.normalize(); 
                refraction = trace(phit - nhit * bias, refrdir, spheres, boxes, depth + 1); 
            } 
            // the result is a mix of reflection and refraction (if the sphere is transparent)
            surfaceColor = ( 
                reflection * fresneleffect + 
                refraction * (1 - fresneleffect) * box->transparency) * box->surfaceColor; 
        } 
        return surfaceColor;
    }
    
} 

void render(const std::vector<Sphere> &spheres, const std::vector<Box> &boxes) 
{ 
    unsigned width = 1920, height = 1080; 
    Vec3f *image = new Vec3f[width * height], *pixel = image; 
    float invWidth = 1 / float(width), invHeight = 1 / float(height); 
    float fov = 45, aspectratio = width / float(height); 
    float angle = tan(M_PI * 0.5 * fov / 180.); 
    // Trace rays
    for (unsigned y = 0; y < height; ++y) { 
        for (unsigned x = 0; x < width; ++x, ++pixel) { 
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio; 
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle; 
            Vec3f raydir(xx, yy, -1); 
            raydir.normalize(); 
            *pixel = trace(Vec3f(0), raydir, spheres, boxes, 0); 
        } 
    } 
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary); 
    ofs << "P6\n" << width << " " << height << "\n255\n"; 
    for (unsigned i = 0; i < width * height; ++i) { 
        ofs << (unsigned char)(std::min(float(1), image[i].x) * 255) << 
               (unsigned char)(std::min(float(1), image[i].y) * 255) << 
               (unsigned char)(std::min(float(1), image[i].z) * 255); 
    } 
    ofs.close(); 
    delete [] image; 
}

int main(int argc, char **argv) 
{ 
    srand48(13); // not entirely sure why we are seeding a random number generator
    std::vector<Sphere> spheres;
    std::vector<Box> boxes; 
    boxes.push_back(Box(Vec3f(2, 2, -10), Vec3f(4,4.0, -15), Vec3f(0.3, 0.3, 0.3), 0, 0));
    // A box cannot be used as a light
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0)); // ground
    spheres.push_back(Sphere(Vec3f( 0.0,      0, -20),     4, Vec3f(1.00, 0.32, 0.36), 1, 0.5)); 
    spheres.push_back(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.80, 0.80, 0.80), 1, 0.0)); 
    spheres.push_back(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.0, 0.77, 0.97), 0, 0.0)); 
    spheres.push_back(Sphere(Vec3f(-5.5,      0, -15),     3, Vec3f(0.90, 0.90, 0.90), 0, 0.0)); 
    // light
    spheres.push_back(Sphere(Vec3f( 0.0,     20, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
    //spheres.push_back(Sphere(Vec3f( 0.0,     10, 10),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(5.0, 3.0, 1.5)));
    render(spheres, boxes); 
 
    return 0; 
} 